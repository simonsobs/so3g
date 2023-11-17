"""HKTree: representation of an HK data archive's structure as a tree
of attributes.

"""

from so3g.hk import getdata
import time
import os
import yaml
import logging


logger = logging.getLogger(__name__)


class HKRef(object):
    """Node in an HKTree attribute tree.

    Because public child attributes are generated dynamically,
    important functionality is hidden in "private members".  The
    ``_load`` and ``_clear`` methods are documented below, but be
    aware also of the following attributes:

    - ``_data``: The loaded data, as a tuple of arrays (t, val).
    - ``_t``: alias for _data[0].
    - ``_val``: alias for _data[1].
    - ``_private``: A dict of information for managing the reference,
      including the full field name, the root tree object, the list of
      child refs.

    """
    _data = None

    def __init__(self, name, parent, terminal):
        super().__init__()
        self._private = {
            'name': name,
            'parent': parent,
            'children': [],
            'terminal': terminal,
            'alias': None,
        }

    @property
    def _t(self):
        if self._data:
            return self._data[0]

    @property
    def _val(self):
        if self.data:
            return self._data[1]

    def __getattr__(self, k):
        # For non-private, invalid attribute names ... return an
        # HKDeadend so user scripts don't immediately fail in cases
        # where a field is not present during this time interval.
        if k[0] == '_':
            raise AttributeError(k)
        return HKDeadend(self._private['name'] + '.' + k,
                         self._private['parent'], False)

    def __repr__(self):
        name = self._private['name']
        post = ' ...' if not self._private['terminal'] else ''
        data = ' [loaded]' if self._data is not None else ''
        return f'<{self.__class__.__name__}:{name}{post}{data}>'

    def _load(self, **kw):
        """Load the data for this node and all child nodes.  The data are
        saved, in each terminal node, in the ``_data`` attribute as a
        (times, values) tuple.

        Returns:
          All the data that was loaded, as a dictionary where the key
          is the full field name and the value is the (times, values)
          tuple.

        """
        return self._private['parent']([self], **kw)

    def _clear(self):
        """Discard any loaded data for this node and all child nodes."""
        if self._private['terminal']:
            self._data = None
        else:
            [x._clear() for x in self._private['children']]


class HKDeadend(HKRef):
    pass


class HKAliasRefs(HKRef):
    def _load(self, **kw):
        return self._private['parent']([self], use_aliases=True, **kw)


class HKTree:
    def __init__(self, start=None, stop=None, config=None,
                 data_dir=None, pre_proc_dir=None,
                 aliases=None, skip=['observatory', 'feeds']):
        """Scan an HK archive, between two times, and create an attribute tree
        representing the HK data.

        Args:

          start (time): Earliest time to include (defaults to 1 day
            ago).
          stop (time): Latest time to include (defaults to 1 day after
            start).
          config (str): Filename of a config file (yaml).  Alternately
            a dict can be passed in directly.
          data_dir (str): The root directory for the HK files.
          pre_proc_dir (str): Directory to use to store/retrieve
            first-pass scanning data (see HKArchiveScanner).
          aliases (dict): Map from alias name to full field name.
            This setting does not override aliases from the config
            file but instead will extend them.
          skip (list of str): Tokens to suppress when turning feed
            names (e.g. "observatory.X.feeds.Y") into tree components
            (e.g. X.Y).

        Notes:
          Initialization of the tree requires a "first pass" scan of
          the HK data archive.  The time range you specify is thus
          very important in limiting the amount of IO activity.  The
          arguments passed are closely related to the load_ranges
          function in this module; see that docstring.

          Config files that work with load_ranges should also work
          here.

        """
        # Parse time args
        now = time.time()
        if start is None:
            start = now - 86400
        else:
            start = getdata.to_timestamp(start)
        if stop is None:
            stop = start + 86400
        else:
            stop = getdata.to_timestamp(stop)

        if aliases is None:
            aliases = {}
        else:
            aliases = dict(aliases)  # copy

        # Use config dict / file?
        if isinstance(config, str):
            config = yaml.safe_load(open(config, 'rb'))
        if config is not None:
            if data_dir is None:
                data_dir = config.get('data_dir')
            if pre_proc_dir is None:
                pre_proc_dir = config.get('pre_proc_dir')

            if config.get('field_list'):
                for k, v in config['field_list'].items():
                    if k not in aliases:
                        aliases[k] = v
            if config.get('skip_tokens'):
                skip = skip + list(config['skip_tokens'])

        # Final default substitutions.
        if data_dir is None:
            data_dir = os.environ['OCS_DATA_DIR']

        # Walk the files -- same approach as load_ranges
        logger.debug('Scanning %s (pre_proc=%s)' % (data_dir, pre_proc_dir))
        hksc = getdata.HKArchiveScanner(pre_proc_dir=pre_proc_dir)
        for folder in range(int(start / 1e5), int(stop / 1e5) + 1):
            base = os.path.join(data_dir, str(folder))
            logger.debug(f' ... checking {base}')
            if not os.path.exists(base):
                continue

            for filename in sorted(os.listdir(base)):
                logger.debug(f' ... ... processing {filename}')
                try:
                    t = int(filename[:-3])
                except ValueError:
                    logger.warning(' ... ... filename does not lead with '
                                   f'timestamp, skipping: {filename}')
                    continue
                if t >= start - 3600 and t <= stop + 3600:
                    hksc.process_file_with_cache(os.path.join(base, filename))

        self._private = {
            'skip': skip,
            'children': [],
        }
        self._private['arc'] = hksc.finalize()
        self._private['fields'] = self._private['arc'].get_fields()

        # Prepare alias reverse map...
        rev_aliases = {v: k for k, v in aliases.items()}
        self._aliases = HKAliasRefs('', self, False)

        # Build attribute tree.
        for k in self._private['fields'][0].keys():
            target = self._find(k, create_missing=True)
            if k in rev_aliases:
                self._add_alias(rev_aliases[k], target)

    def _find(self, name, create_missing=False):
        """Find the HKRef in the attribute tree associated with a full field
        name.  If create_missing is False, this will return None if
        the name is not found in the tree; if create_missing is True,
        any missing HKRef will be created and the terminal node
        returned.

        """
        tokens = name.split('.')
        target = self
        name = ''
        for i, t in enumerate(tokens):
            end_point = (i == len(tokens) - 1)
            name += t
            if t not in self._private['skip']:
                t = t.replace('-', '_')
                if not hasattr(target, t) or \
                   isinstance(getattr(target, t), HKDeadend):
                    if not create_missing:
                        return None
                    new_node = HKRef(name, self, terminal=end_point)
                    target._private['children'].append(new_node)
                    setattr(target, t, new_node)
                target = getattr(target, t)
            name += '.'
        return target

    def _add_alias(self, alias, full_ref, _strip_prefix=None):
        """Add a field to the set of aliases.

        Args:
          alias (str): The alias key.
          full_ref (str or HKRef): The full name of the field or an
            HKRef for the field.

        Notes:

          If the full_ref is a non-terminal HKRef, then all children will
          be added in.  The alias for each will be constructed by
          combining the provided alias string and the attribute names
          of each child node, joined with '_'.

        """
        if isinstance(full_ref, str):
            target = self._find(full_ref)
            if target is None:
                raise ValueError(f'No node found for "{full_ref}"')
        else:
            target = full_ref

        # Handle non-terminal nodes
        name = target._private['name']
        if not target._private['terminal']:
            if _strip_prefix is None:
                _strip_prefix = name
            for t in target._private['children']:
                self._add_alias(alias, t, _strip_prefix=_strip_prefix)
            return
        if target._private['alias'] is not None:
            logger.warning('Field {k} is already aliased as {k}, '
                           'ignoring new alias.')
            return
        if _strip_prefix:
            assert(name.startswith(_strip_prefix))
            alias = alias + '_' + name[len(_strip_prefix)+1:]
        target._private['alias'] = alias
        clean_alias = alias.replace('.', '_').replace('-', '_')
        setattr(self._aliases, clean_alias, target)
        self._aliases._private['children'].append(target)

    def __getattr__(self, k):
        if k[0] == '_':
            raise AttributeError(k)
        return HKDeadend(k, self, False)

    def _load(self, **kw):
        """Call _load() on all child attributes and return a giant data dictionary
        (see HKRef._load).

        """
        return self(self._private['children'], **kw)

    def _clear(self):
        """Drop all loaded data."""
        for c in self._private['children']:
            c._clear()

    def __call__(self, fields, use_aliases=False, **kw):
        def get_targets(fs):
            if isinstance(fs, HKRef):
                if fs._private['terminal']:
                    return [fs]
                return get_targets(fs._private['children'])
            fields = []
            for f in fs:
                fields.extend(get_targets(f))
            return fields
        field_targets = get_targets(fields)
        # Tuples (field_name, key_name)
        names = [f._private['name'] for f in field_targets]
        if use_aliases:
            names = [(n, n if f._private['alias'] is None
                      else f._private['alias'])
                     for n, f in zip(names, field_targets)]
        else:
            names = [(n, n) for n in names]
        data = self._private['arc'].simple([n[0] for n in names])
        # Populate the endpoints.
        for t, vects in zip(field_targets, data):
            t._data = vects
        return {k[1]: d for k, d in zip(names, data)}
