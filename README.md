# SubgraphContainment
## Introduction
A subgraph query finds all data graphs in a graph
database each of which contains the given query graph. Existing
work takes the indexing-filtering-verification (IFV) approach
to first index all data graphs, then filter out some of them
based on the index, and finally test subgraph isomorphism on
each of the remaining data graphs. This final test of subgraph
isomorphism is a sub-problem of subgraph matching, which finds
all subgraph isomorphisms from a query graph to a data graph.
As such, in this paper, we study whether, and if so, how to
utilize efficient subgraph matching to improve subgraph query
processing. Specifically, we modify leading subgraph matching
algorithms and integrate them with top-performing subgraph
query processing algorithms. Our experimental results show
that (1) the slow verification method in existing IFV algorithms
can lead us to over-estimate the gain of filtering; and (2) The
latest subgraph matching algorithms can be easily modified
to the vertex connectivity based filtering-verification subgraph
query processing algorithms, which is competitive with the IFV
algorithms while do not rely on any indices. Considering the
problems in the indexing of IFV algorithms such as the scalability
issue and the update cost with data graphs modification, the
second finding is very important as it provides an alternative
approach of the IFV paradigm, which can scale up subgraph
queries to hundreds of thousands of data graphs and graphs of
tens of thousands of vertices, and makes an initial step towards
the indexing-free subgraph query processing.

For the details, please refer to our ICDE'2019 paper (to be appeared)
"Scaling Up Subgraph Query Processing with Efficient Subgraph Matching"
by [Shixuan Sun](https://github.com/shixuansun) and [Dr. Qiong Luo](http://www.cse.ust.hk/~luo/).
If you have any further questions, please feel free to contact us.

Please cite our paper, if you use our source code.

* "Shixuan Sun and Qiong Luo. Scaling Up Subgraph Query Processing with
   Efficient Subgraph Matching. ICDE 2019."



## Compile
Under the src directory, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Execute
Execute the binary with the following command ./SubgraphQuery -d data_graphs -df FB -q query_graphs -qf FB,
in which -d specifies the input of the data graphs and -q specifies the input of the query graphs.
The -df and -qf configure the input format of data graphs and query graphs respectively.
Please do not change them.

Example:

```zsh
./SubgraphQuery -d ../../dataset/test_data_graphs.graph -df FB -q ../../dataset/test_query_graphs.graph -qf FB
```

## Input
Both the input data sets and query sets can be a collection of graphs.
Each graph starts with 't # ID' where ID is the name of the graph. A vertex and an edge are formatted
as 'v VertexID LabelId' and 'e VertexId VertexId 0' respectively. Note that we require that the vertex
id is started from 0 and the range is [0,|V| - 1] where V is the vertex set. The following
is an input sample. You can also find a sample data sets and query sets under the dataset folder.

Example:

```zsh
t # 0
v 0 0
v 1 3
v 2 3
v 3 0
v 4 0
e 0 1 0
e 0 2 0
e 0 3 0
e 3 4 0
t # 1
v 0 0
v 1 0
v 2 4
v 3 3
v 4 3
e 0 1 0
e 0 2 0
e 1 3 0
e 1 4 0
```

## Experiment Datasets
The real world datasets and the corresponding query sets used in our paper can be downloaded [here](https://hkustconnect-my.sharepoint.com/:u:/g/personal/ssunah_connect_ust_hk/Ed2ElrYUVV9AkVjXm9otZe8BajhfU5N6G3mKRb7rv0kXBw?e=HyfhiF).
As the synthetic datasets are large, we do not upload them. You can easily generate them by following the instruction in our paper.

## Acknowledgement
We would like to thank Prof. Xuemin Lin, Dr. Fei Bi (the authors of [1]), Dr. Nikos Ntarmos, Dr. Foteini Katsarou (the authors of [2]), Dr. James Cheng (the author of [3]),
Dr. Karsten Klein, Dr. Nils Kriege (the authors of [4]), Dr. Rosalba Giugno and Dr. Dennis Shasha (the authors of [5]) for kindly sharing their
source code and answering our questions.

* [1]. F. Bi, L. Chang, X. Lin, L. Qin, and W. Zhang. Efficient subgraph
matching by postponing cartesian products. In SIGMOD, 2016.
* [2]. F. Katsarou, N. Ntarmos, and P. Triantafillou. Performance and scala-
bility of indexed subgraph query processing methods. In PVLDB, 2015.
* [3]. J. Cheng, Y. Ke, W. Ng, and A. Lu. Fg-index: towards verification-free
query processing on graph databases. In SIGMOD, 2007.
* [4]. K. Klein, N. Kriege, and P. Mutzel. Ct-index: Fingerprint-based graph
indexing combining cycles and trees. In ICDE, 2011.
* [5]. R. Giugno, V. Bonnici, N. Bombieri, A. Pulvirenti, A. Ferro, and
D. Shasha. Grapes: A software for parallel searching on biological
graphs targeting multi-core architectures. In PloS one, 2013.