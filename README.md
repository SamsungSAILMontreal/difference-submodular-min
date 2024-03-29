# Difference of Submodular Minimization via DC Programming


Code to reproduce results of the paper [Difference of Submodular Minimization via DC Programming](https://arxiv.org/pdf/2305.11046.pdf)

## To reproduce results in the paper
- run dsm_job.sh script 
- plot results using dsm_plot.m

## Datasets
- Datasets are already included in datasets folder
- The speech dataset was included with the code from [1].
- The mushroom dataset was downloaded from https://www.openml.org/search?type=data&status=active&id=24. The train/test splits  used in our experiments were generated by running the python script mushroom.py (requires NumPy, SciPy, and pandas).

```bash
cd datasets 
python mushroom.py
```
# Citation
```
@InProceedings{elhalabi2023dsm,
      title={Difference of Submodular Minimization via DC Programming}, 
      author={Marwa El Halabi and George Orfanides and Tim Hoheisel},
      booktitle = {Proceedings of the 40th International Conference on Machine Learning},
      year={2023},
}
```
## Acknowledgements
- We use some functions from [1] (mainly the implementation of MNP algorithm).
- We use for plotting a modified version of the padcat function from [2].
- We use for plotting the distinguishable_colors function from [3].

[1]  Francis Bach, Matlab Submodular package (version 2.0), https://www.di.ens.fr/~fbach/submodular/. Retrieved September, 2016.

[2] Jos (10584) (2023). PADCAT (https://www.mathworks.com/matlabcentral/fileexchange/22909-padcat), MATLAB Central File Exchange. Retrieved October 30, 2022.

[3] Tim Holy (2023). Generate maximally perceptually-distinct colors (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors), MATLAB Central File Exchange. Retrieved May 17, 2021.
