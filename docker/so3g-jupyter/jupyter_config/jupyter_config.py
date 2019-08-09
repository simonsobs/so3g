import os

c.NotebookApp.ip = '0.0.0.0'
c.NotebookApp.port = 8888

c.NotebookApp.allow_root = True
c.NotebookApp.open_browser = False

config_dir = os.path.dirname(__file__)
pw_file = os.path.join(config_dir, '.jupyter_pw_hash')
with open(pw_file, 'r') as f:
    c.NotebookApp.password = f.read().strip()

c.NotebookApp.password_required = True

c.NotebookApp.notebook_dir = u'/config/'

