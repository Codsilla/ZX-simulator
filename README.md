# ZX-simulator
Classical simulator of quantum computing with stabilser decompositions through ZX-calculus.
This is a re-implementation of the simulator of [quizx](https://github.com/Quantomatic/quizx) but now including graphs cuts (using KaHyPar hypergraph partionning libray).

This is a quick and dirty implementation as it was developed as I was experimenting.

As of now the fastest decomposition function is `decomp_most_connected_cat` with `imbalence` set as 0.5.

Note that you'll need to have a working build of KahyPar with bindings for rust which you can find [here](https://github.com/tuomas56/kahypar-rs). You will then have to update the dependencies to include your path to KaHyPar-rs.
