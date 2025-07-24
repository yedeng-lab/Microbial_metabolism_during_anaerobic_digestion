#bipartite network analysis
#use sparcc
fastspar --iterations 20 --exclude_iterations 10 --threshold 0.3 --otu_table ./otutable.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
mkdir ./bootstrap_counts
fastspar_bootstrap --otu_table ./otutable.txt --number 100 --prefix ./bootstrap_counts/fake_data
mkdir bootstrap_correlation
for j in `seq 0 99`
do
fastspar --iterations 20 --exclude_iterations 10 --threshold 0.3 --otu_table ./bootstrap_counts/fake_data_$j.tsv --correlation bootstrap_correlation/cor_fake_data_$j.tsv --covariance bootstrap_correlation/cov_fake_data_$j.tsv
done
fastspar_pvalues --otu_table ./otutable.txt --correlation median_correlation.tsv --prefix ./bootstrap_correlation/cor_fake_data_ --permutations 100 --outfile pvalues.tsv
done

