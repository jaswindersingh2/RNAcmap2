# RNAcmap2
An improved fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis


|![](../docs/figure_1.png)
|----|
| <p align="center"> <b>Figure 1:</b> The architecture of the RNAcmap2 pipeline. CSS: Consensus Secondary Structure. CM: Covariance Model. L: Length of the input RNA sequence.|

1. `conda activate venv_rnacmap2`

2. `./run.py --neff no_hit --dca_method gremlin`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.000 	 0.000 		 0.000		    19		  0.0    
direct_infernal     0.000 	 0.000 		 0.000		    19		  0.0    
RNAcmap             0.000 	 0.000 		 0.000		    19		  0.0    
RNAcmap_meta        0.112 	 0.145 		 0.097		    19		  1.0    
RNAcmap2_meta       0.132 	 0.180 		 0.110		    19		  1.4
```

3. `./run.py --neff low --dca_method gremlin`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.006 	 0.007 		 0.006		    85		  0.0    
direct_infernal     0.153 	 0.196 		 0.127		    85		  2.0    
RNAcmap             0.169 	 0.218 		 0.140		    85		  2.1    
RNAcmap_meta        0.214 	 0.274 		 0.181		    85		  4.1    
RNAcmap2_meta       0.386 	 0.499 		 0.322		    85		  13.0
```

4. `./run.py --neff median --dca_method gremlin`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.036 	 0.047 		 0.029		    31		  1.0    
direct_infernal     0.357 	 0.486 		 0.284		    31		  17.0    
RNAcmap             0.411 	 0.554 		 0.334		    31		  26.5    
RNAcmap_meta        0.461 	 0.617 		 0.383		    31		  31.6    
RNAcmap2_meta       0.532 	 0.708 		 0.436		    31		  96.1
```

5. `./run.py --neff high --dca_method gremlin`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.074 	 0.104 		 0.058		    114		  2.2    
direct_infernal     0.589 	 0.805 		 0.467		    114		  321.4    
RNAcmap             0.633 	 0.862 		 0.502		    114		  636.5    
RNAcmap_meta        0.630 	 0.857 		 0.501		    114		  605.1    
RNAcmap2_meta       0.630 	 0.858 		 0.501		    114		  605.1
```

6. `./run.py --neff all --dca_method gremlin`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.040 	 0.056 		 0.032		    249		  1.0    
direct_infernal     0.366 	 0.496 		 0.293		    249		  15.6    
RNAcmap             0.398 	 0.538 		 0.319		    249		  27.5    
RNAcmap_meta        0.427 	 0.574 		 0.346		    249		  35.3    
RNAcmap2_meta       0.497 	 0.665 		 0.402		    249		  100.6
```

7. `./run.py --neff all --dca_method plmc`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.008 	 0.011 		 0.007		    249		  1.0    
direct_infernal     0.000 	 0.000 		 0.000		    249		  0.0    
RNAcmap             0.398 	 0.538 		 0.320		    249		  27.5    
RNAcmap_meta        0.436 	 0.586 		 0.354		    249		  35.3    
RNAcmap2_meta       0.506 	 0.676 		 0.412		    249		  100.6
```

8. `./run.py --neff all --dca_method mfdca`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.054 	 0.074 		 0.043		    249		  1.0    
direct_infernal     0.000 	 0.000 		 0.000		    249		  0.0    
RNAcmap             0.416 	 0.559 		 0.337		    249		  27.5    
RNAcmap_meta        0.452 	 0.606 		 0.368		    249		  35.3    
RNAcmap2_meta       0.525 	 0.700 		 0.427		    249		  100.6
```

9. `./run.py --neff all --dca_method plmdca`


```
		       F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.025 	 0.034 		 0.020		    249		  1.0    
direct_infernal     0.000 	 0.000 		 0.000		    249		  0.0    
RNAcmap             0.408 	 0.549 		 0.330		    249		  27.5    
RNAcmap_meta        0.444 	 0.595 		 0.362		    249		  35.3    
RNAcmap2_meta       0.518 	 0.691 		 0.423		    249		  100.6
```
