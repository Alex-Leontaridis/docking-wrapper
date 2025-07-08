from config import Config
config = Config()
config.mgltools_path = get_path("/test/path")
config.output_dir = get_path("/test/output")
import json
with open('test_config.json', 'w') as f:
    json.dump(config.get_config_dict(), f)
print('Configuration saved')

with open('test_config.json', 'r') as f:
    loaded = json.load(f)
print(f"MGLTools path: {loaded['mgltools_path']}")
print(f"Output dir: {loaded['output_dir']}") 