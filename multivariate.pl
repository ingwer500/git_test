#!/usr/bin/perl
#===============================================================================
##
##         FILE: multivariaten.pl
##
##        USAGE: ./multivariaten.pl --help
##
##  DESCRIPTION: statistical calculations extra
##
##      OPTIONS: see help
## REQUIREMENTS: Statistics::Descriptive from  Schlomi Fish, FileHandle
##         BUGS: no plausibility tests; if the deviation equals to null ( by same distribution) the script will be interrupted.
##        NOTES: do not enter after the last line in data file please ! there is no plausibilty test. Look at the source code
##       AUTHOR: ingwer ( ingwer@freeshell.org ),
## ORGANIZATION: 
##      VERSION: 0.03
##      CREATED: 22.07.2013 14:31:54
##     REVISION: alpha
##===============================================================================
use strict ;
use warnings ;
use POSIX ;
use Statistics::Descriptive ;
use FileHandle ;

use Getopt::Long ;
#-------------
my $interpolation = 'linear' ;

my $optfile ;
my @optbeta ;
my @optgamma ;
my @optinterval ;
my @optdescriptive ;
my @optsxyvalue ;
my @optquantile ;
my @optcovariance ;
my @optmxyvalue ;
my @optpearson ;
my @optcorrelation ;
my @optitest ;
my @optchii ;
my @optinterpol ;
my @opthelp ;
my @optverbose ;
#----------------
my @matrix ;
my $rows ;
my $cols ;

my $result = GetOptions(
						"file=s"							=> \$optfile,
						"Sxy|S"								=> \@optsxyvalue ,
						"covariance|c"						=> \@optcovariance,
						"mvalue|m"							=> \@optmxyvalue,
						"pearson|p"							=> \@optpearson,
						"correlation|r"						=> \@optcorrelation,
						"itest"								=> \@optitest,
						"chiInd"							=> \@optchii,
						"quantiles|q"						=> \@optquantile ,
						"descriptive|d"						=> \@optdescriptive,
						"intervall=s{2}"					=> \@optinterval, 
						"gamma=s{3}"						=> \@optgamma,
						"verbose"							=> \@optverbose,
						"help|h"							=> \@opthelp,
				);

			
sub __Help__{
	print "Usage: $0 [ --file file ] 
			 [--intervall value1 value2 : aproximated ]
			 [--gamma 0|value1 inf|value2 df : aproximated ]
			 [-S|--Sxy]
			 [-p|--pearson]
			 [-c|--covariance]
			 [-m|--mvalue]
			 [-r|--correlation]
			 [-q|--quantiles]
			 [--descriptive]
			 [--chiInd: Chi for Indepedence test]
			 [--itest : Indepedence test]
			 [--verbose]\n" ;
	exit ;
}

#-----
my $fh = FileHandle->new() ;

sub __FileOpen__ {

	my @data ;
	my $file = shift @_  ;

	if (  $fh->open("< $file") ) {
		for ( <$fh> ) {
			push (@data, [ split( /\s/ , $_ ) ] ) ;		# make ( /\s/,$_ ) it better please !
		}
		$fh->close ;
	}
return @data ;

}


	if ( $optfile ){
		@matrix = __FileOpen__( $optfile );
	}else{ __Help__ ; };

		$rows=@matrix ;			#
		$cols=@{$matrix[0]} ;	#


#------
#
#	my @matrix=(
#		[2450736.89448,138192781.124, 56108274.9375, -389.28934253, -5.5997843, 19.015439, 0.65538606],
#		[2450736.90559, 138187717.479, 56128219.4051, 397.539983932, -5.2958679, 22.037272, 0.88950438],
#		[2450736.9167, 138182444.931, 56150022.0529, 1254.04754107, -5.6898297, 23.228155, 0.88555444],
#		[2450736.93893, 138170907.338, 56195689.4192, 2903.07201826, -6.2759552, 24.179035, 0.83261846],
#		[2450736.98337, 138145523.29, 56289970.7907, 5950.35778361, -6.8659225, 24.803789, 0.76229525],
#		[2450737.02781, 138118525.721, 56385731.7476, 8796.96438698, -7.16868, 25.042354, 0.72361985],
#		[2450737.1167, 138062062.567, 56578955.6999, 14175.7631663, -7.4983983, 25.243222, 0.68211734],
#		[2450737.29448, 137944082.466, 56967883.0369, 24330.4470959, -7.8254248, 25.369863, 0.64561973],
#		[2450737.47226, 137822297.315, 57357874.9963, 34102.0826381, -8.0205678, 25.40276, 0.62844001],
#		[2450737.82781, 137571446.223, 58138266.4374, 53112.0806181, -8.2960737, 25.394701, 0.61147173],
#		...) ;
#	my $rows=@matrix ;			#
#	my $cols=@{$matrix[0]} ;	#
#
#------


my $stat=Statistics::Descriptive::Full->new() ;
my $statexp=Statistics::Descriptive::Full->new() ;

sub z_s($$) {

	my $interpolation = shift ;
	my $i = shift ;
	$stat->clear() ;
	$statexp->clear() ;

	my ( $mean, $count, $sum, $var, $min, $max, $deviation, $summ2, @fit) ;
	
	
	my $j = 0 ;
	my @extract ;
	$summ2 = 0; 
	for my $rows ( @matrix ) {
	if (  ! defined @$rows[ $i ] ) { printf "%s%d\n","Ungleiche Anzahl der Proben. Moegliche Zeile ", $j+1; exit  ; }

		elsif ( $interpolation eq 'linear' ){	# 	
			my $xi = @$rows[$i] * @$rows[$i] ;
			$summ2 = $summ2 + $xi ;				# sum Xi**2 linear
			$extract[$j++] =  @$rows[$i] ;		#
			$stat->add_data( @$rows[$i] ) ;

			}
			#elsif( $interpolation eq 'logarithmic' ){ do something for logarithmic } # test
	}
				( $mean, $count, $sum, $var, $min, $max, $deviation, @fit) = ( 
				$stat->mean(), 
				$stat->count, 
				$stat->sum(),							# not used 
				$stat->variance, 
				$stat->min, 							# not used
				$stat->max,								# not used
				$stat->standard_deviation(),
				$stat->least_squares_fit); 				# not used
				my @samplerange = $stat->sample_range ;	# not used
				my @quantilen = ( $stat->quantile(0),$stat->quantile(1),$stat->quantile(2),$stat->quantile(3),$stat->quantile(4) ) ;

				my @frequency;
				push ( @frequency,$stat->frequency_distribution( \@quantilen ) ) ;	# not used

				my $var2 = $var * ( $count - 1 ) / ( $count ) ;				#
				my $deviation2 = sqrt($var2) ;								# 
				return $mean , $count , $sum , $var , $var2, $deviation, $deviation2, $summ2 ,\@extract, \@fit, \@samplerange, \@quantilen, \@frequency ;

}


my @colvalue = () ;

sub __colvalue__($) {
	my $interpolation = shift ;
	for ( 0..$cols-1 ) { push ( @colvalue , [ z_s ( $interpolation, $_ ) ] ) ; };

}
__colvalue__($interpolation) ;

sub __colvaluechi__($) {
	my $interpolation = shift ;
	my @colvalueforchi=() ;
	for ( 0..$cols-1 ) { 
		push ( @colvalueforchi , [ z_s ( $interpolation, $_ ) ] ) ; 
	};
	return @colvalueforchi ; 
}




#---- start one column case ---

	if ( $cols == 1 ) {
		print qq/\x79\x6f\x75\x20\x63\x61\x6e\x20\x77\x72\x69\x74\x65\x20\x69\x74\x20\x73\x65\x6c\x66\x0a/ ;
		#start point print "@{$colvalue[0]}","\n" ;
	   	exit ;	
			
	}

#---- end one column case ---

sub __column__ {
	my @column;
	my $i=0 ;
	my @colvalue = @colvalue ;
	for my $col (  @colvalue ) {
		$column[$i] = $col->[8] ;
		$i++ ;
	} ;
	return @column ;
}

sub xixn($$$) {
	my @column = __column__ ;
	my $a=0 ;
	my $interpol = shift ;
	my $b=shift ;
	my $j=shift ;
	my $x1xn=0 ;

		for ( @matrix ) {
			if( $interpol eq 'linear'){	
				my $var = $column[$j][$a] * @{$matrix[$a]}[$b] ;
			 	$x1xn = $var + $x1xn ;
				$a++ ;
			}
		} ;
		 return $x1xn ;					# guter Wert

} ;



sub __xixn__($) {				# function , which delivers $i und $j values into the xixn() and returns the array @x1xn.
	my @x1xn ;				#
	my $interpol= shift ;
	for ( my $i = 1 ; $i <= $cols-1 ; $i++  ){
			for ( my $j = $i ; $j <= $cols-1 ; $j++ ) {
			 push ( @x1xn,  xixn( "$interpol",$j,$i-1)  ) ;
		}	

	}
	return @x1xn ;
}
#__xixn__( $interpolation );

sub  __S_xixn__ {			# F(Xi, Xk) = Summ[ Xi*Xk ] - n * mean[ Xi ] * mean[ Xk ]  values.

my @colvalue= @colvalue ;
my @S_xixn	= ();
my $var		= 0 ;
my $a		= 0 ;
my @x1xn  	= __xixn__($interpolation)		;	# Hier ist der Einsprungspunkt fuer Interpolierung
my $n = $colvalue[0][1] ;
for ( my $i = 1 ; $i <= $cols-1; $i++) {
		for( my $j = $i ; $j <= $cols-1 ; $j++ ) {

				$var = $x1xn[ $a++ ] -  $n * $colvalue[$i-1][0] * $colvalue[$j][0] ;
				push ( @S_xixn, $var  ) ;
		}
	}
	return @S_xixn ;
}

#---ChiSquare Start---
#
sub __zeile__{
	my @zeilen ;
	for ( @matrix ) {
			push ( @zeilen, $_ ) ;
	}
	return @zeilen ;

}

my @zeilensummen;

sub __zeilensumme__ {
	my @zlist = __zeile__ ;
	#my @zsummen;
	for ( my $i=0 ; $i < $rows ; $i++ ) { 
		my $zsumm=0  ;
		for ( my $a=0 ; $a < $cols ; $a++ ) {
			my $var = $zlist[$i][$a] ;
			$zsumm = $zsumm + $var ;
			}
			push (@zeilensummen , $zsumm ) ;
		}
		#return  @zsummen ;
}

my @spaltensummen ;

sub __spaltensumme__ {
	my @werte = __colvaluechi__('linear')  ;
	my @swerte = @werte  ;
	for ( @swerte ) {

		push(  @spaltensummen , $$_[2]  );

	}
	#return @ssummen ;
}


my $summederSsummen=0 ;
sub __summederSsummen__ {
	my $summederSsummen=0 ;
			for ( @spaltensummen ) {
	
	 		$summederSsummen = $summederSsummen + $_ ;
		}
		#return $summederSsummen;
}

my $summederZsummen=0 ;
sub __summederZsummen__ {
		for ( @zeilensummen ) {
	
			 $summederZsummen = $summederZsummen + $_ ;
		}

		#return $summederZsummen ;
}

my @produktZsummenSsummen=() ;
sub __prodZsumSsum__ {
		for my $ssumm ( @spaltensummen ){

			for my $zsumm ( @zeilensummen ) {
				my $produktZsummenSsummen = $zsumm * $ssumm ;
				push ( @produktZsummenSsummen, $produktZsummenSsummen );
			}
		}
		#return @produktZsummenSsummen ;
}

sub __prodZsumSsumdurchN__ {		# 

	my @prodZsumSsumdurchN ;

		for ( @produktZsummenSsummen ) {
			push ( @prodZsumSsumdurchN , $_ / $summederZsummen  ) ;
		   	
		}
	return @prodZsumSsumdurchN ;
}

#---- Start Unabhaengigkeitstest ( Kontingenzanalyse / Kreuztabelle)  ----
__zeilensumme__;__spaltensumme__;__summederSsummen__;__summederZsummen__;__prodZsumSsum__ ;
sub __ChiQuadratWerte__ {



	my @zlist=__zeile__ ;
	my @prodZsumSsumdurchN =  __prodZsumSsumdurchN__  ;
	my $Freiheitsgrade = ($rows-1) * ($cols-1) ;
	my $l=0 ;
	my @zuerwartendeMatrix ;
	my @ChiQuadratWerte ;

	for ( my $i=0 ; $i < $cols; $i++ ) {
			for(  my $j=0 ; $j < $rows;  $j++ ) {
				if ( $Freiheitsgrade > 1 ) {

	if ( $prodZsumSsumdurchN[$l] == 0 ) {  print "no Chi" ; return undef ; } ;
	  my $chi = ( ( $zlist[$j][$i] - $prodZsumSsumdurchN[$l] ) * ( $zlist[$j][$i] - $prodZsumSsumdurchN[$l] ) ) / $prodZsumSsumdurchN[$l] ;
	  $ChiQuadratWerte[$l] = $chi ;
	  $zuerwartendeMatrix[$i][$j] = $prodZsumSsumdurchN[$l] ;
	  $l++ ;
	 }else{
	if ( $prodZsumSsumdurchN[$l] == 0 ) {  print "no Chi" ; return undef ; } ;
	my $chi = ( ( abs( $zlist[$j][$i] - $prodZsumSsumdurchN[$l]) - 0.5 ) * ( abs( $zlist[$j][$i] - $prodZsumSsumdurchN[$l]) - 0.5 ) ) / $prodZsumSsumdurchN[$l] ;
	 $ChiQuadratWerte[$l] = $chi ;
	 $zuerwartendeMatrix[$i][$j] = $prodZsumSsumdurchN[$l] ;
	 $l++ ; }

	 } ;

	}
	return \@ChiQuadratWerte, \@zuerwartendeMatrix ;
}


sub __PrintzuerwartendeMatrix__ {

	my @chiquadratwerte =__ChiQuadratWerte__ ;
	my @zuerwartendeMatrix = @{$chiquadratwerte[1]} ;
	my $Freiheitsgrade = ($rows-1) * ($cols-1) ;
	print "Zu erwartende Matrix Freiheitsgrade $Freiheitsgrade\n";
	print "-"x50,"\n" ;
		for ( my $i=0; $i<$rows ; $i++) {
			for ( my $j=0; $j<$cols ; $j++ ) {
			printf "%.4f ","$zuerwartendeMatrix[$j][$i] ";
			}
		print "\n";
			}
		print "-"x50,"\n";
}


sub __ChiQuadrat__ {
		my @ChiQuadratWerte = __ChiQuadratWerte__ ;
		my $ChiSum=0 ;
			for ( @{$ChiQuadratWerte[0]} ) {
			 $ChiSum = $ChiSum + $_ ;
		}
		return $ChiSum ;
}

sub __Covariance_XiXk__ {
		my @Covariance;
		my @colvalue = @colvalue ;
		my $n = $colvalue[0][1] ;
			for my $Sxixn ( __S_xixn__ ) {
			push(@Covariance , $Sxixn / $n ) ;

	}
	return @Covariance ;
}

sub __Print_SummeXiXk__ {

	my $a=0 ;
	my @SummeXiXk = __S_xixn__ ;

	for ( my $i = 0 ; $i < $cols  ; $i++ ) {
			for ( my $j = $i+1 ; $j < $cols ; $j++ ) {
				print  "S[ x$i,x$j ] =  $SummeXiXk[ $a++ ]", "\n";
				}
	}

}


sub __PrintCovariancen__(@) {
	my $a=0 ;
	my @colvalue =  @_ ;
	my $n = $colvalue[0][1] ;
	my @Covariancen =  __Covariance_XiXk__ ;
	for  (  my $i = 0 ; $i<$cols ; $i++  ) {
		for  ( my $j=$i+1   ; $j<$cols ; $j++ ) {
			print "Cov( x$i,x$j ) = $Covariancen[ $a ]", "\n" ;
			$a++ ;
				}
		}
}


sub __Mxywerte__(@) {
	
	my @Mxywerte =() ;
	my @Sxy = __S_xixn__;
	my @Sxx = @_ ;
	my $n = $Sxx[0][1] ;
	my $a = 0 ;
	my $k = 1;
	for( my $i = 0 ;  $i < $cols ; $i++ ){
		for( my $j = $i ; $j < $cols-1 ; $j++ ) {
				my $m1 = ( $Sxy[$a] ) / ( $Sxx[ $i ][ 4 ] * $n  ) ;	
				my $m2 = ( $Sxy[$a] ) / ( $Sxx[ $j+1 ][ 4 ] * $n )  ;
				my $var1 = sprintf "%s%d%s%d%s%d%s%d%s%0.8f","m[ ",$i,".",$j+1,",",$i,".",$i," ] = ",$m1 ;
				my $var2 = sprintf "%s%d%s%d%s%d%s%d%s%0.8f","m[ ",$i,".",$j+1,",",$j+1,".",$j+1," ] = ", $m2 ;
				push ( @Mxywerte ,join ("\t\t", $var1, $var2 ) ) ;
				$a++ ;
			}
		};

	return @Mxywerte ;
};


# beide sind kopien von einander.
sub __KorrelationsWerteNachPearson__(@) {
	
	my @Sxy = __S_xixn__;
	my @Sxx = @_ ;
	my $n = $Sxx[0][1] ;
	my $a = 0 ;
	for( my $i = 0 ;  $i < $cols ; $i++ ){
		for( my $j = $i ; $j < $cols-1 ; $j++ ) {
				my $m1 = ( $Sxy[$a] ) / ( $Sxx[ $i ][ 4 ] * $n ) ;
				my $m2 = ( $Sxy[$a] ) / ( $Sxx[ $j+1][ 4 ] * $n)  ;
				my $rxixkquadrat = $m1 * $m2  ;
				print "rpearson [ x$i,"," x",$j+1," ] = ",sqrt ( $rxixkquadrat ),"\n" ;
				$a++ ;
			}
		}
	}
#---- stop ----

sub __PrintMxyWerte__ {
		my $n=__column__;
		my $anzahlderMwerte=( 0.5 * ($n-1) * ($n-1) + 0.5 * ($n-1)  ) ;
		my @array=__Mxywerte__(@colvalue) ;
		 for( my $i = 0 ; $i < $anzahlderMwerte  ; $i++ ) {
			print "$array[$i]","\n";
		}
}


sub __PrintKorrelationen__(@) {

	my @covariancen = __Covariance_XiXk__ ;
	my @colvalue = @_ ;
	my $a = 0 ;
	for ( my $i = 0     ;    $i < $cols     ; $i ++ ) {

			if(  $colvalue[$i][4]  == 0 ){ print "Mindestens Probe $i ist gleichverteilt.\n" ; exit ; }

		for ( my $j=$i+1 ;    $j < $cols     ; $j++) {

			if( $colvalue[$i][4]  == 0 || $colvalue[$j][4] == 0 ){ print "Mindestens Probe $j ist gleichverteilt.\n" ; exit ; }
			print "r[ x$i,x$j ] = ", $covariancen[ $a ] / sqrt( $colvalue[$i][4] * $colvalue[$j][4]), "\n" ;
			$a++ ;

		}
	}
}


sub __Quantilen__(@) {

	for my $fit ( @_ ) {
		print "Q[0%] < $$fit[11]->[0], Q[25%] < $$fit[11]->[1], Q[50%] < $$fit[11]->[2], Q[75%] < $$fit[11]->[3], Q[100%] < $$fit[11]->[4]",	"\n" ;
	}
}

#-- Start Standard-Normal Verteilung --
sub __GaussVerteilung__(@) {

	my $h = 0.0001 ;	# is a good value
	my $pi = 4*atan(1) ;
	
	my ($a,$b) = @_  ;

	if( ! ($a =~ /^([0-9]+(\.[0-9]+)?)$/ ||	$b =~ /^([0-9]+(\.[0-9]+)?)$/ ) ){	__Help__ ;	} ;
	if ( $a <= -10 || $b >=10  ) { print "Value between -10 10 please" ; exit ; } ;

	my $n = abs ( $b-$a ) / $h ;
	my $coef = ( $b - $a ) / ( 3 * $n ) ;
	my ( $sum,$sum4,$sum2 ) = ( 0,0,0 ) ;
	my $fa = ( 1 / sqrt( 2 * $pi) ) * exp( -0.5 * ( $a * $a ) ) ;
	my $fb = ( 1 / sqrt( 2 * $pi) ) * exp( -0.5 * ( $b * $b ) ) ;
			for my $i ( 1..$n ) {

			my $zi = $a + $h * ( $i - 1) ;
				if ( $i % 2 != 0 ) {
					my $fxi = 4 * ( ( 1 / sqrt( 2 * $pi) ) * exp( -0.5 * ( $zi ) * ( $zi ) ) ) ;
					$sum4 = $sum4 + $fxi ; }
				else{
					my $fxj = 2 * ( ( 1 / sqrt( 2 * $pi) ) * exp( -0.5 * ( $zi) * ( $zi ) ) ) ;
					$sum2 = $sum2 + $fxj ;
				}
			}
			$sum = $coef * ( $fa + $sum4 + $sum2 + $fb ) ;
			my $result=sprintf "%.4f", $sum ;
			return $result ;
}


sub __Gamma__($$$) {

	# If you increase $h=0.1 you get better "n!"- values ( at example ./multivariate --file data.dat --gamma 0 inf 100 delivers better result )
	# but worse sqrt( PI ) value ( at example ./multivariate --file data.dat --gamma 0 inf 0.5 delivers worse result )
	# i decided for better sqrt(PI) value 
	my $h = 0.0001 ;	# is a good value for sqrt(PI)
	#my $h = 0.1 ;		# is a good value for n! 
	#--my $pi = 4*atan(1) ;

	my ($a,$b,$v) = @_  ;		# von 0 bis inf mit Freiheitsgrad v

		#if ( $a  == 0  ) { $a = 0.0001 ;	}	# simulate Zerro good value.
		if ( $a  == 0  ) { $a = 0.00005 ;	}	# simulate Zerro better value !
		if ( $b eq 'inf' ) { $b = 180   ;	}	# simulate Infinite
		if (  !$a || !$b || !$v  ) {
			__Help__ ;
		}
		
	# here comes Taylor aproximation
	my $n = abs ( $b-$a ) / $h ;
	my $coef = ( $b - $a ) / ( 3 * $n ) ;
	my ( $sum,$sum4,$sum2 ) = ( 0,0,0 ) ;
	my $fa = ( ( exp( -$a ) ) * $a ** ( $v - 1  ) ) ;
	my $fb = ( ( exp( -$b ) ) * $b ** ( $v - 1  ) ) ;
			for my $i ( 1..$n ) {

			my $zi = $a + $h * ( $i - 1) ;
				if ( $i % 2 != 0 ) {
					my $fxi =  4 *  ( ( exp(  -$zi ) ) * ( $zi **  ( ( $v - 1 ) ) ) )  ;
					$sum4 = $sum4 + $fxi ; }		#
				else{
					my $fxj =  2 * ( ( exp(  -$zi ) ) * ( $zi **  ( ( $v - 1 ) ) ) ) ;
					$sum2 = $sum2 + $fxj ;
				}
			}
			$sum = $coef * ( $fa + $sum4 + $sum2 + $fb ) ;
			my $result=sprintf "%.4f", $sum ;		# optimal 4f
			return  $result ;
}

#------Gamma2 Stop------

sub __ChiVerteilung__($$$) {

	# Chi Indepedence Test is Not a homogeniety test  !
	my $h = 0.0001 ;	# is a good value 
	my $pi = 4*atan(1) ;
	my ($a,$b,$v) = @_  ;		# von 0 bis inf mit Freiheitsgrad v
	my $gamma = __Gamma__( 0, 50 , $v * 0.5 ) ;	# 50 guter Wert je hoeher n!- Wert desto genauer, kompromiss zw. Tabelle und Taschenrechner

	my $n = abs ( $b-$a ) / $h ;
	my $coef = ( $b - $a ) / ( 3 * $n ) ;
	my ( $sum,$sum4,$sum2 ) = ( 0,0,0 ) ;
	my $fa = ( 1 / ( 2 ** ( $v * 0.5 ) * $gamma ) ) * ( ( exp( -$a * 0.5 ) ) * $a ** ( ( $v * 0.5 ) - 1  ) ) ;
	my $fb = ( 1 / ( 2 ** ( $v * 0.5 ) * $gamma ) ) * ( ( exp( -$b * 0.5 ) ) * $b ** ( ( $v * 0.5 ) - 1  ) ) ;
			for my $i ( 1..$n ) {

			my $zi = $a + $h * ( $i - 1) ;
				if ( $i % 2 != 0 ) {
					my $fxi =  4 * (  ( 1 / ( 2 ** ( $v * 0.5 ) * $gamma) ) * ( ( exp(  -0.5 * $zi ) ) * $zi **  ( ( $v*0.5 ) - 1 ) ) ) ;
					$sum4 = $sum4 + $fxi ; }		#
				else{
					my $fxj =  2 * (  ( 1 / ( 2 ** ( $v * 0.5 ) * $gamma) ) * ( ( exp(  -0.5 * $zi ) ) * $zi **  ( ( $v*0.5 ) - 1 ) ) ) ;
					$sum2 = $sum2 + $fxj ;
				}
			}
			$sum = $coef * ( $fa + $sum4 + $sum2 + $fb ) ;
			my $result=sprintf "%.4f", $sum ;		# 4 decimal places
			return  $result ;
}
my $Chi = __ChiQuadrat__ ;
my $v = ( $rows - 1 ) * ( $cols - 1 )  ;	#  Freiheitsgrade des Unabhaengigkeitstests

#- Stop Chi Verteilung
sub __KonfidenzIntervall__(@) {

	my ($a,$b) = @_ ;
	my $i=0;
	my @data=@colvalue ;	#  
	my $n = $data[0][1] ;
	my $posib = __GaussVerteilung__( $a,$b ) ;
		for my $var ( @data ) {
		 my $deviation = $var->[5] ;
		 my $sterror = $deviation / sqrt ( $n ) ;
		 printf "%s%d%s%.4f%s%.4f%s%.4f%s\n","p ( x",$i," = ",$posib," ) = [ " ,$var->[0] - $posib*$sterror," <= mu <=  ",$var->[0] + $posib*$sterror," ]" ;
		 $i = $i + 1 ;
			}
}

sub __PrintDescriptive__(@) {
	my $i=0 ;
	my @data=@colvalue ;
	my @descriptive = () ;
	for my $var ( @data ) { 

		my $var1 = sprintf "%s%d %s %.4f   %s%d %s %.4f   %s%d\n","mean[ x",$i,"] = ",$var->[0],"s[ x",$i,"] = " ,$var->[5],"n = ",$var->[1];
		print "$var1" ;
		$i=$i+1 ;
		}
}

	
	if ( @optgamma ){
		my ( $value1, $value2, $df ) = @optgamma ;
		print "Gamma ( ","a = ", $value1 ,";","b = ",$value2,";","df = ",$df ," ) = ",__Gamma__( $value1,$value2,$df ),"\n"  ; }
	if( @optsxyvalue ){
		__Print_SummeXiXk__ ;
	}
	if( @optcovariance ){
		__PrintCovariancen__( @colvalue ) ;
	}
	if(  @optmxyvalue){
		__PrintMxyWerte__ ;
	}
	if( @optinterval ){
		__KonfidenzIntervall__( @optinterval ) ;
	}
	if( @optquantile ){
		__Quantilen__( @colvalue ) ;
	}
	if( @optpearson ){
		__KorrelationsWerteNachPearson__( @colvalue );
	}
	if( @optcorrelation ){
		__PrintKorrelationen__( @colvalue );
	}
	if( @optdescriptive ){
		__PrintDescriptive__(@colvalue) ;
	}
	if( @opthelp ){
		__Help__ ;
	}
	if( @optverbose ){
		print "interpolation: $interpolation\n" ;
	}
	if( @optchii ){
		 printf "%s%.4f\n","Chisquare = ", __ChiQuadrat__ ;
	}
	if( @optitest ){
		printf "%s%.4f%s%d%s%.4f\n","f ( chi = ", $Chi," v = ", $v,") = ",__ChiVerteilung__(0, $Chi , $v ),"\n" if ( $Chi <= 90 )  ; 
		# indepedence test makes sense by chi Value smaller than 90
	}

#---------------------------------- Stop Hauptteil---------------------------------------#
