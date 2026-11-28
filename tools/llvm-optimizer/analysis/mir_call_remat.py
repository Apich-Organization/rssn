"""Find post-call broadcasts and reloads needed on only one Win64 machine edge.

Input is LLVM X86 MIR after virtregrewriter. Static matches are candidates for
inspection. The script does not rewrite MIR or model branch execution frequency.
"""

import argparse
import collections
import re
from pathlib import Path


BLOCK = re.compile(r"^  bb\.(\d+)(?:[. (][^:\n]*)?:\s*$", re.M)
SUCCESSOR = re.compile(r"%bb\.(\d+)")
CALL = re.compile(r"\bCALL64pcrel32 [@&]([^,\s]+), csr_win64\b")
BROADCAST = re.compile(
    r"^(?:renamable )?\$([yz]mm\d+) = V(?:P)?BROADCAST\w+ "
    r"\$rip, 1, \$noreg, (%const\.\d+), \$noreg\b.*from constant-pool"
)
CALL_LIKE = re.compile(r"^(?:CALL|TAILJMP|TCRETURN)")
STACK_RESTORE = re.compile(
    r"^(?:renamable )?\$(\w+) = (?:MOV64rm|VMOV\w+rm) "
    r"(%(?:fixed-)?stack\.\d+), 1, \$noreg, 0, \$noreg\b.*\bload "
)
FRAME_OBJECT = re.compile(
    r"^  - \{ id: (\d+), .*?type: (\S+),[^}]*}", re.M
)
FIXED_LOAD = re.compile(r"\bload .* from (%fixed-stack\.\d+)\b")
SPILL_LOAD = re.compile(r"\bload .* from (%stack\.\d+)\b")


def machine_functions(text):
    for document in re.split(r"^---\s*$", text, flags=re.M):
        name = re.search(r"^name:\s*(\S+)\s*$", document, re.M)
        body = re.search(r"^body:\s*\|\s*\n", document, re.M)
        if name and body and re.search(r"^noVRegs:\s+true$", document, re.M):
            fixed = re.search(r"^fixedStack:\s*\n(.*?)(?=^\w)", document, re.M | re.S)
            local = re.search(r"^stack:\s*\n(.*?)(?=^\w)", document, re.M | re.S)
            slots = {}
            if fixed:
                for slot in FRAME_OBJECT.finditer(fixed[1]):
                    if "isImmutable: true" in slot[0]:
                        slots[f"%fixed-stack.{slot[1]}"] = "fixed-stack reads"
            if local:
                for slot in FRAME_OBJECT.finditer(local[1]):
                    if slot[2] == "spill-slot":
                        slots[f"%stack.{slot[1]}"] = "spill reloads"
            yield name[1], document[body.end():], slots


def machine_blocks(body):
    headers = list(BLOCK.finditer(body))
    blocks = {}
    for index, header in enumerate(headers):
        end = headers[index + 1].start() if index + 1 < len(headers) else len(body)
        successors, liveins, code = [], set(), []
        for line in body[header.end():end].splitlines():
            if line.startswith("    successors:"):
                successors = SUCCESSOR.findall(line)
            elif line.startswith("    liveins:"):
                liveins = set(re.findall(r"\$(\w+)", line))
            elif line.startswith("    "):
                code.append(line.strip())
        blocks[header[1]] = (successors, liveins, code)
    return blocks


def register_aliases(reg):
    if re.fullmatch(r"[yz]mm[0-9]+", reg):
        return {prefix + reg[1:] for prefix in "xyz"}
    if re.fullmatch(r"r(?:8|9|1[0-5])", reg):
        return {reg + suffix for suffix in ("", "d", "w", "b")}
    families = (
        ("rax", "eax", "ax", "ah", "al"), ("rcx", "ecx", "cx", "ch", "cl"),
        ("rdx", "edx", "dx", "dh", "dl"), ("rbx", "ebx", "bx", "bh", "bl"),
        ("rsi", "esi", "si", "sil"), ("rdi", "edi", "di", "dil"),
        ("rbp", "ebp", "bp", "bpl"),
    )
    return next((set(family) for family in families if reg == family[0]), set())


def scan(blocks, slots, callee):
    counts = collections.Counter()
    matches = []
    for block, (successors, _, code) in blocks.items():
        for index, line in enumerate(code):
            call = CALL.match(line)
            if not call or (callee and callee != call[1]):
                continue
            counts["direct calls"] += 1
            end = next((pos for pos in range(index + 1, len(code))
                        if CALL_LIKE.match(code[pos])), len(code))
            window = code[index + 1:end]
            counts["post-call constant-pool broadcasts"] += sum(
                bool(BROADCAST.match(inst)) for inst in window)
            counts["post-call fixed-stack reads"] += sum(
                bool(FIXED_LOAD.search(inst)) for inst in window)
            counts["post-call local stack reads"] += sum(
                bool(SPILL_LOAD.search(inst)) for inst in window)
            if end != len(code) or len(successors) != 2:
                continue
            for definition, inst in enumerate(window):
                if broadcast := BROADCAST.match(inst):
                    reg, source = broadcast.groups()
                    kind = "broadcasts"
                elif (reload := STACK_RESTORE.match(inst)) and reload[2] in slots:
                    reg, source = reload.groups()
                    kind = slots[source]
                else:
                    continue
                aliases = register_aliases(reg)
                if not aliases:
                    continue
                mention = re.compile(r"\$(?:" + "|".join(sorted(aliases)) + r")\b")
                if any(mention.search(tail) or source in tail
                       for tail in window[definition + 1:]):
                    continue
                for advance, exit_to in (successors, successors[::-1]):
                    _, advance_live, advance_code = blocks[advance]
                    _, exit_live, _ = blocks[exit_to]
                    if advance_live & aliases or reg not in exit_live:
                        continue
                    for following in advance_code:
                        if mention.search(following):
                            break
                        if CALL_LIKE.match(following):
                            if (next_call := CALL.match(following)) and next_call[1] == call[1]:
                                counts[f"edge-only {kind}"] += 1
                                matches.append((block, advance, exit_to, call[1],
                                                reg, source))
                            break
    return counts, matches


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mir", type=Path)
    parser.add_argument("--function", required=True, help="machine function name substring")
    parser.add_argument("--callee", help="exact symbol of the Win64 direct call")
    args = parser.parse_args()
    found = False
    for name, body, slots in machine_functions(args.mir.read_text()):
        if args.function not in name:
            continue
        found = True
        counts, matches = scan(machine_blocks(body), slots, args.callee)
        print(name)
        for block, advance, exit_to, callee, reg, source in matches:
            print(f"  bb.{block}: bb.{advance} -> @{callee}, "
                  f"bb.{exit_to} needs ${reg} from {source}")
        for label, count in counts.items():
            print(f"  {label}: {count}")
    if not found:
        parser.error(f"no machine function matching {args.function!r}")


if __name__ == "__main__":
    main()
