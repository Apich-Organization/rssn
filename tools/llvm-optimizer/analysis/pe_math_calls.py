"""Locate candidates for imported libm calls in an x86-64 PE release image.

Byte matches are leads. Confirm instruction boundaries in a full function
disassembly before using a candidate as evidence.
"""

import argparse
import bisect
import collections
import re
import struct
from pathlib import Path

MATH = frozenset(("pow", "exp", "exp2", "expm1", "log", "log1p", "cosh", "sinh"))


class PE:
    def __init__(self, path):
        self.data = path.read_bytes()
        header, = self.read("<I", 0x3C)
        if self.data[header:header + 6] != b"PE\0\0\x64\x86":
            raise ValueError("expected an x86-64 PE image")
        sections, = self.read("<H", header + 6)
        option_size, = self.read("<H", header + 20)
        optional = header + 24
        if self.read("<H", optional) != (0x20B,):
            raise ValueError("expected a PE32+ optional header")
        self.base, = self.read("<Q", optional + 24)
        self.directories = [self.read("<II", optional + 112 + i * 8) for i in (1, 3)]
        self.sections = {}
        for i in range(sections):
            entry = optional + option_size + 40 * i
            name = self.data[entry:entry + 8].split(b"\0", 1)[0]
            virtual, rva, size, start = self.read("<IIII", entry + 8)
            self.sections[name] = (rva, start, min(virtual, size))

    def read(self, fmt, offset):
        return struct.unpack_from(fmt, self.data, offset)

    def offset(self, rva):
        for base, start, size in self.sections.values():
            if base <= rva < base + size:
                return start + rva - base
        raise ValueError(f"unmapped RVA 0x{rva:x}")

    def imports(self, symbols=MATH):
        result = {}
        rva, length = self.directories[0]
        if not length:
            return result
        start = self.offset(rva)
        for p in range(start, start + length, 20):
            lookup, _, _, _, address = self.read("<IIIII", p)
            if not address:
                break
            table = self.offset(lookup or address)
            for index in range((len(self.data) - table) // 8):
                item, = self.read("<Q", table + 8 * index)
                if not item:
                    break
                if item >> 63:
                    continue
                offset = self.offset(item) + 2
                name = self.data[offset:self.data.index(b"\0", offset)].decode("ascii")
                if name in symbols:
                    result[address + index * 8] = name
        return result

    def boundaries(self):
        rva, length = self.directories[1]
        if not length:
            return []
        start = self.offset(rva)
        return sorted(self.read("<III", start + i)[:2] for i in range(0, length - 11, 12))

    def calls(self, symbols=MATH):
        imports = self.imports(symbols)
        base, start, size = self.sections[b".text"]
        end = start + size
        thunks = {}
        hits = []
        for match in re.finditer(rb"\xff[\x15\x25]", self.data[start:end]):
            at = start + match.start()
            if at + 6 > end:
                continue
            rva = base + at - start
            delta, = self.read("<i", at + 2)
            symbol = imports.get(rva + 6 + delta)
            if symbol is None:
                continue
            if self.data[at + 1] == 0x25:
                thunks[rva] = symbol
            else:
                hits.append((rva, symbol))
        if thunks:
            for match in re.finditer(rb"\xe8", self.data[start:end]):
                at = start + match.start()
                if at + 5 > end:
                    continue
                rva = base + at - start
                delta, = self.read("<i", at + 1)
                if symbol := thunks.get(rva + 5 + delta):
                    hits.append((rva, symbol))
        return sorted(hits)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dll", type=Path)
    parser.add_argument("--limit", type=int, default=12)
    parser.add_argument("--symbol", action="append", help="restrict imports to these names")
    args = parser.parse_args()
    pe = PE(args.dll)
    boundaries = pe.boundaries()
    starts = [start for start, _ in boundaries]
    regions = collections.defaultdict(list)
    calls = pe.calls(set(args.symbol) if args.symbol else MATH)
    for rva, symbol in calls:
        index = bisect.bisect_right(starts, rva) - 1
        if index >= 0 and rva < boundaries[index][1]:
            regions[boundaries[index]].append((rva, symbol))
    print(f"{len(boundaries)} unwind regions, {len(calls)} candidate call sites")
    for (start, end), sites in sorted(
        regions.items(),
        key=lambda item: (len({s for _, s in item[1]}), len(item[1])),
        reverse=True,
    )[:args.limit]:
        counts = collections.Counter(symbol for _, symbol in sites)
        print(
            f"RVA 0x{start:x}..0x{end:x} (VA 0x{pe.base + start:x}..0x{pe.base + end:x})",
            dict(sorted(counts.items())),
        )
        print("  " + " ".join(f"{symbol}@0x{rva:x}" for rva, symbol in sites))


if __name__ == "__main__":
    main()
