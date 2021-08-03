# NormRcode
R Code for Normalization Of Methylation Array with Meffil based on Dr. Drew Prat's workflow and master metadata for biowulf environment. 
Takes metadata and extra samples (only centrix ids in separate file) requares normalization otherwise will extract prebuild betavalues and run
<br>

## Examples
To run pipeline in noninteractive mode on large memory nodes. 
* Run "Placenta" samples from master database output results to: /data/MDATA/TRANSFER/NewBATCHdir and re-normalize placenta samples with group of extra samples listed in: TRANSFER/newbatch.txt. That text file have one sample per line first column only requared with id like this: 9934625026_R06C02
"Placenta:NewCNSdir:TRANSFER/newbatch.txt"

~~~
sbatch /data/MDATA/NormRcode/norm.sh Placenta:NewBATCHdir:TRANSFER/newbatch.txt Placenta:NewBATCHdir:TRANSFER/newbatch.txt
~~~

* Run "Placenta" samples from master database output results to: /data/MDATA/TRANSFER/NewBATCHdir no normalization just reuse beta values stored. Output umap file. 
~~~
sbatch /data/MDATA/NormRcode/norm.sh Placenta:NewBATCHdir Placenta:NewCNSdir:
~~~

* Run "Placenta" and "Gastrointestinal" samples together from master database output results to: /data/MDATA/TRANSFER/NewBATCHdir. Renormalize everything and Output umap file for R-Shiny. 
~~~
sbatch /data/MDATA/NormRcode/norm.sh Placenta|Gastrointestinal:NewBATCHdir:TRANSFER/newbatch.txt 
~~~

* Run ALL samples together from master database output results to: /data/MDATA/TRANSFER/NewBATCHdir. Renormalize everything and Output umap file for R-Shiny. You can have a single sample in  TRANSFER/newbatch.txtto make it happened.

~~~
sbatch /data/MDATA/NormRcode/norm.sh .:NewBATCHdir:TRANSFER/newbatch.txt 
~~~


