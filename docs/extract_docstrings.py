import os
import sys
import inspect
import re


def install_fake_spt3g(dest_dir):
    dest_dir = os.path.join(dest_dir, 'spt3g/core')
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    with open('%s/__init__.py' % dest_dir, 'w') as fout:
        for core_class in [
                'G3Frame',
                'quat',
                'G3VectorQuat',
                'G3VectorDouble',
                'G3VectorInt',
                'G3VectorString',
                'G3VectorBool',
        ]:
            fout.write('class %s:\n    pass\n' % core_class)


def install_fake_so3g(dest_dir, src_dir=None):
    if src_dir is None:
        src_dir = '%s/python' % dest_dir
    dest_dir = os.path.join(dest_dir, 'so3g')
    if not os.path.exists(dest_dir):
        os.symlink(src_dir, dest_dir)


def create_docstring_shells(module, fout, select=None):
    def print(text='', end='\n'):
        fout.write(text + end)

    print('# This file is generated automatically by scanning a compiled\n'
          '# C++ boost-python module, to facilitate documentation builds.\n'
          '# It should not be present in the master branch!\n'
          '# Edits will be lost.\n')

    for cname, cls in module.__dict__.items():

        if select is not None and cname not in select:
            continue

        if cname == 'version':
            print('def version():')
            print('  return "%s"' % cls())

        elif inspect.isclass(cls):
            print('class %s:' % cname)
            if hasattr(cls, '__doc__'):
                print('  """%s"""' % cls.__doc__)
            print('  pass')
            for k, v in cls.__dict__.items():
                if k[0] == '_':
                    continue

                if not inspect.isroutine(v):
                    continue
                # Determine staticness from the __dict__ version.
                is_static = isinstance(v, staticmethod)

                # Get the docstring from the getattr version.  (The
                # __dict__ version doesn't work for staticmethods.)
                docstring = getattr(cls, k).__doc__

                # Attempt to extract the call signature from the first
                # non-empty line of the docstring.
                doclines = docstring.strip('\n').split('\n')
                sig_line = ''
                if len(doclines):
                    sig_line = doclines[0]
                m = re.fullmatch(
                    r'(\w+)\( ' + r'(.*)' + r'\) -> (.*) :', sig_line)
                if m is None:
                    arglist = ['*args', '**kwargs']
                else:
                    _, arglist, ret_type = m.groups()
                    arg_toks = re.findall(r'\((\w+)\)(\w+)', arglist)
                    arglist, optional = [], False
                    for t, n in arg_toks:
                        if t == '':
                            optional = True
                        elif optional:
                            arglist.append('%s=something' % n)
                        else:
                            arglist.append(n)

                    # Rewrite docstring.
                    arg_ofs = 1
                    if is_static:
                        arg_ofs = 0
                    doclines[0] = ('%s(' % k
                                   + ', '.join(arglist[arg_ofs:])
                                   + ') -> %s' % ret_type)

                if is_static:
                    print('  @staticmethod')
                print('  def %s():' % (k))
                print('    """' + '\n'.join(doclines) + '"""')
                print('    pass')
            print()
        else:
            print('# No handler for "%s"\n' % cname)


def prepare_readthedocs(src_branch='master',
                        dest_branch='readthedocs',
                        dest_file=None,
                        version_file=None):
    import so3g
    if dest_file is None:
        docs_dir = os.path.split(__file__)[0]
        dest_file = os.path.join(
            docs_dir, '../python/_libso3g_docstring_shells.py')

    if version_file is None:
        version_file = os.path.join('../docs/_so3g_rtd_version.txt')

    for cmd in [
            'git checkout %s' % dest_branch,
            'git pull --ff-only',
            'git merge --no-ff %s' % src_branch,
            None,
            'git add %s' % dest_file,
            'git add %s' % version_file,
            'git commit --allow-empty -m "Doc-ready build : %s"' % so3g.version(),
    ]:
        if cmd is None:
            # Write the lib shell and version file.
            print('Creating docstring shells in %s...' % dest_file)
            create_docstring_shells(so3g, open(dest_file, 'w'))
            print('Writing version to %s' % version_file)
            open(version_file, 'w').write(so3g.__version__)
        else:
            print('run: %s' % cmd)
            code = os.system(cmd)
            if code != 0:
                print(' FAILED with code %i' % code)
                return 10
    print()
    print('Merging into the %s branch seems to have succeeded. Check that \n'
          'everything is right and then push to github.  You are on branch:\n'
          % dest_branch)
    os.system('git symbolic-ref HEAD --short')
    print()
    return 0


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--select', action='append', default=None)
    parser.add_argument('-o', '--output-file', default=None)
    parser.add_argument('-p', '--prep-rtd', action='store_true')
    parser.add_argument('--source-branch', default='master')
    parser.add_argument('--dest-branch', default='readthedocs')
    args = parser.parse_args()

    if args.prep_rtd:
        sys.exit(prepare_readthedocs(
            args.source_branch, args.dest_branch, args.output_file))

    else:
        import so3g
        if args.output_file is not None:
            fout = open(args.output_file, 'w')
        else:
            fout = sys.stdout
        create_docstring_shells(so3g, fout, select=args.select)
