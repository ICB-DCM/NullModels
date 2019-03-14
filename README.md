# NullModels

The aim of a null models is to predict the expected effect of two drugs *A* and *B* given as the combined dose pair *(a, b)*. The single response curves are here modeled as hill curves.

The following models are supported:

- [x] Loewe
- [x] Bliss
- [x] Hand
- [x] Tallarida (given as lower and upper bound due to non uniqueness)
- [x] Highest Single Agent

## Organization of the code:

- **main.m**: Run this script to read in data + evaluate data as well as to replicate the figures from [our publication](https://www.nature.com/articles/s41598-019-38907-x). You need to download the data from the supplementary of [O'Neil et al. 2016](http://mct.aacrjournals.org/content/15/6/1155.long) and change the values of `singleDrugData =...` (line 6) and `CombinationData = ...`(line 9) to the directories of the corresponding files on your coumputer.

### Data structure
The data structure is given in the corresponding folder.

- **Drug**: All functionality corresponding Hill Curves (+ their fitting) and the single response data.
- **Combination**: Contains two drugs and gives all functionality  corresponding to the null models and stores the combination treatment data.
- **Cell Line**: Contains Drug + Combination for one particular cell line (Hill Curves are fitted for each cell line separately).
- **Data**: Stores all data, that is read in (data comes from [O'Neil et al. 2016](http://mct.aacrjournals.org/content/15/6/1155.long)),

### Plots
The functions are in the folder **Plots**
- **IsobolePlot**: Isoboles of the different null models. Used in Figure 7. (Call `IsobolePlot(D, 2, 3) ` ).
- **CorrPlots.m**: Correlates the predicion of different null models. Used in Figure 8. (Call `CorrPlots(D.CellLines{3})`) .
- **VolumeMetricConceptPlot** and **VolPlot** Plot volumes between the predicted and the measured response surfaces. (Not used in the publication. Call e.g. `VolumeMetricConceptPlot` and `VolPlot(D)`).

## Publication

[Sinzger, Vanhoefer, Loos, Hasenauer (2019) Comparison of null models for combination drug therapy reveals Hand model as biochemically most plausible. Scientific Reports, 9, 3002](https://www.nature.com/articles/s41598-019-38907-x)

If you have any questions regarding the implementations and/or need help: We encourage you to contact us! 
