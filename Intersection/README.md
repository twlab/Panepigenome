## Graph-based interval-overlap analysis

This workflow identifies overlaps between CpG and variant intervals directly in pangenome graph coordinates.

CpG and variant intervals are represented as oriented graph paths with terminal offsets. Variants are indexed by the graph segments they traverse, allowing each CpG interval to be compared only with variants sharing at least one segment. Overlap is then calculated at the segment-base level while accounting for graph orientation.

Reported results include the CpG identifier, variant identifier, overlap length and source interval annotations.

This approach enables interval intersections to be evaluated directly in graph space without conversion to a single linear reference coordinate system.
