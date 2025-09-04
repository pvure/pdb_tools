import json
import re

def jsonc_to_json(jsonc_file, json_file):
    with open(jsonc_file, 'r') as f:
        content = f.read()
    
    # Remove single-line comments
    content = re.sub(r'//.*', '', content)
    
    # Remove multi-line comments
    content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
    
    # Parse and reformat JSON
    data = json.loads(content)
    
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)

jsonc_to_json('/Users/pranayvure/Downloads/binder_partial_sample.jsonc', '/Users/pranayvure/Downloads/binder_partial_sample.json')