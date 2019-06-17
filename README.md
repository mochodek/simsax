# Similarity-based on Symbolic Aggregate approXimation (SimSAX)
SimSAX: A Measure of Project Similarity Based on Symbolic Approximation Method and Software Defect Inflow

Cite:
```
@article{ochodek2019simsax,
  title={{SimSAX}: A Measure of Project Similarity Based on Symbolic Approximation Method and Software Defect Inflow},
  author={Ochodek, Miros{\l}aw and Staron, Miroslaw and Meding, Wilhelm},
  journal={Information and Software Technology},
  doi={10.1016/j.infsof.2019.06.003}
}
```

## Installation

You can use pip to install the package once it is cloned:

```
git clone https://github.com/mochodek/simsax.git
cd simsax
pip install -e .
```

You can go through the provided Jupyter Notebooks to learn how to use the SimSAX calculator and run simulation of the meausure parameters n, w, a. To run the calculator you need to have two CSV files containing time series for two projects (see the data-examples folder) and run the simsax tool. The tool takes the following parameters:
* A directory where the input files are stored
* A csv input filename with a backlog for the first project
* A csv input filename with a backlog for the second project
* A column name for the time series to be considered
* Window length
* Word length
* Alphabet size

Example:
```
simsax "./data-examples/" "eclipse_jdt_sample_backlog.csv" "eclipse_platform_sample_backlog.csv" "inflow_all" 64 5 6
```
