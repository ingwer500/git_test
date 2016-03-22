Perl script for statistical calculations
========================================

- [multivariate.pl](http://ingwer500.github.io "ingwer on github.com")

1. Produce file with numbers
	* Example
[ user@box ]$ for i in 1 2 3 4 5 6 7 8 9 ; do echo $RANDOM $RANDOM $RANDOM $RANDOM  done > data.dat

2. Run the script
	* Example
[ user@box ]$ perl ./multivariate.pl --file data.dat --descriptive --intervall -1.96 1.96

3. You can See help by typping
[ user@box ]$ perl ./multivariate --help
