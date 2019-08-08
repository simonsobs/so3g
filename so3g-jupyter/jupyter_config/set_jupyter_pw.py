from notebook.auth import passwd
import os

pw_hash = passwd()
config_dir = os.path.join(os.environ["OCS_CONFIG_DIR"], "jupyter")
pw_file = os.path.join(config_dir, '.jupyter_pw_hash')
with open(pw_file, 'w+') as f:
    f.write(pw_hash)
