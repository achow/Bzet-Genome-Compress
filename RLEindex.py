from RLE import rle
import sys

indexes = sys.stdin.readlines()
bitstring = rle(None)

for bit in indexes:
    bitstring.SET(int(bit))

print(bitstring.size())
