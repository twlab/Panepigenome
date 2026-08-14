## Projection of individual assembly coordinates to graph coordinates

Individual assembly coordinates were projected to HPRC2 graph coordinates using haplotype-specific walks encoded in the GFA representation of the Minigraph-Cactus graph, with GRCh38 used as the backbone.

For each chromosome, the corresponding haplotype walk was parsed into an ordered series of oriented graph segments, and cumulative segment offsets were calculated along the walk. Input intervals, including CpG intervals, were then mapped to this cumulative coordinate system.

Each projected interval was represented by:

- the traversed graph path;
- the start offset within the first graph segment;
- the end offset within the last graph segment; and
- the total interval length.

For intervals spanning multiple graph segments, the complete intervening graph path was retained. Segment orientation was explicitly considered, and sequences traversing reverse-oriented segments were reverse-complemented when reconstructing interval sequences.

Projected interval lengths and sequences were compared with the corresponding intervals in the source assemblies to verify coordinate consistency.

### Input

- HPRC2 Minigraph-Cactus GFA files
- Haplotype-specific assembly intervals, including CpG intervals

### Output

- Graph paths corresponding to each input interval
- Start and end offsets within graph segments
- Projected interval lengths
- Reconstructed interval sequences for validation

## Contact

For questions about the lcpg2graph.py script, please contact:

**Wenjin Zhang**
Email: [wenjin@wustl.edu](mailto:wenjin@wustl.edu)
