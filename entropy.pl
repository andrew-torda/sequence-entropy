#!/usr/bin/perl
# 16 March 2003
# 4 nov 2014. This version does a very primitive calculation of the probability
# of one sequence within an alignment if the -f flag is given.
# at each site.
# There is a bug in the -m option (domains). Gaston triggered an undefined
# variable warning at line "    $$domain_changes[$len - 1] = 0;" in do_domains.
# Take the output from a blast / psi-blast search.
# Extract some data for plotting from the alignments.
# rcsid = $Id: entropy.pl,v 1.32 2018/06/26 09:44:36 torda Exp $
=pod

=head1 NAME

entropy.pl

=head1 SYNOPSIS

B<entropy.pl> S<[ B<-e> I<num>]> S<[ B<-o> I<outfile> ] >
S<[B<-a>  S<B<-s> I<seq_name>> I<attfile> [B<-t> I<5 | 6 | 20 | 21>]]> S<[ B<-s> I<num>]> I<alignment>

=head1 DESCRIPTION

Eat blast output in default format.  Collect some data for multiple
sequence alignments.

=head1 OPTIONS

=over

=item B<-a> I<attfile>

Write a file of output in the format for a chimera attribute file. By
default, use 20 symbols and ignore gaps. See B<-t> below.

=item B<-d> I<squashed_seq_file>

Often we are interested in an alignment from the point of view of one
sequence. If we look at this sequence in the middle of an alignment it
is full of gaps. The B<-d> option will take each sequence and print
out only the positions which are present in the query sequence as
defined in the B<-f> option. This option only works if the B<-f>
option has been used. The squashed output will be written to
I<squashed_seq_file>.

=item B<-e> I<number>

Alignments will only be included if their B<e_value> is less than
I<number>.  A default value of 0.005 is used, following the psi-blast
defaults.

=item B<-f> I<name>

The input is fasta format. All the sequences, with
gaps, are the same length. I<name> is the name of the protein you are
interested in. The intention is that you have made a fancy multiple
sequence alignment with something like maffta or t_coffee.

The string I<name> must be present in exactly one of
the sequences. This means that your file might contain a line like

 >gi|194338707|gb|ACF49281.1| CSZ2 [Eimeria acervulina]

and you would say C<-f CSZ2>. I<name>, however, must be unique. If it
is not, the program will print out all occurrences and stop. If you
have duplicates, you probably want to remove them. Otherwise, use a
more unique string in the argument.

If you use the B<-f> flag, it is assumed that you are not working with
blast, so other flags are ignored. For example , the B<-e> flag has no
effect.

As of Nov 2014, this version will calculate the probability of a
sequence within the context of the members of the alignment. It is
primitive. See the note below.


=item B<-m>

Imagine you are looking for doMains (note the capital M) within a
multiple sequence alignment. Some pieces of sequence are present in a
subset of sequences. You want to find these and call them
domains. If you turn on this option, we will look at each site in the
alignment and see in which sequences it is present. We count the
number of changes between positions.
More formally, we visit each column in the alignment and make a bit vector saying if a site is present or absent. We print out the hamming distance between sucessive sites.
The output goes into the plot file as yet another column,
so there is not argument to this option.

When we see a change between columns I<i> and I<i+1>, we add to the
total for both I<i> and I<i+1>. That is, domain 1 runs up to I<i>. The
next domain starts at I<i+1>.

=item B<-o> I<outplotfile>

Data for plotting will be written to outplotfile.
Default is stdout.

=item B<-s> I<offset>

Residue numbering in PDB files often does not start with one, or the
files do not agree with sequence entries. Use this option to correct
for this.  The value I<offset> will be added to all residue numbers on
output.  For example, your PDB file starts with residue 16. Your
sequence starts from the 16th residue and a program like I<chimera>
labels residues starting from 16. Then say <-s 15> to all 15 to all
output.

=item B<-t> I<5|6|20|21>

Only used with the B<-a> option. If you are going to print out an
attribute file for chimera, you chould use

=over

=item 20
Normal entropy with 20 amino acids and ignore gaps.

=item 21
Treat gaps as a 21st amino acid

=item 5
Use the alphabet reduced to five amino acid classes.

=item 6
Use the reduced alphabet, but treat gaps as a sixth residue.

=back

=back

=head1 OUTPUT

We write to a file (or stdout) with a file for plotting with

  res       number        entropy  entropy   frction gaps [residue]
  number    homologues    20 sym   21 sym    present num  [name]

=over

=item S<res number>

Is the residue number in the alignment.

=item S<number homologues>

Is the number of homologues found at this position in the alignment.

=item S<entropy 20 sym>

Is the entropy, but with a 20 symbol alphabet. It does not include gaps.

=item S<entropy 21 sym>

Is the sequence entropy at this position, calculated using a 21
symbol alphabet. It treats gaps as a character.

=item S<gaps num>

The number of gaps within aligned sequences at this position.

=item S<gaps fraction>

The number of gaps at this position, divided by the number of sequences
aligned at this position.

=back

=head1 NOTES AND OPERATION

The script happily reads the output from a multi-round psiblast
calculation. It walks down the file, looking for the start of
each round. At the end, it calls C<fpos()> to go back to the
start of the last search.

If the query aligns to more than one part of the "subject" sequence,
we will only see the first.  After detecting an alignment, the code
looks for the next line beginning with a I<<> character which is the
sign of the start of the next alignment.

Usage is now more often on alignments in fasta format.
This is usually going to be better than just using the blast output.
This leads to a reasonable workflow:

=over

=item * run blast search and set output to just sequence namess

=item * use blast's "blastdbcmd" to get the full sequences that you want

=item * generate a careful alignment with a program like mafft

=item * run this script on the alignment and probably use the S<B<-f abcde>>
option to tell the script that everything should be relative to
your original sequence, B<abcde>.

=back

The output is now normalised by using the logarithm, base the number
of symbols.

=head2 Implementation notes

If the B<-f> flag is set, we are dealing with fasta format and have a specific sequence of interest. Usually this ends up in $q_seq and $seq_ndx will tell us that we have a special sequence. If, we have fasta format, but no particular sequence, we will still have $q_seq, but the variable $seq_ndx should be undef(). This tells later parts of the code not to print out things like residue compatibility.

=head1 TO DO

=over

=item * When reading fasta format, add a check to see if the input
sequences with gaps have the correct length.

=item * Currently, leading gaps in a fasta file seem to get marked as
C<not_char> rather than gaps. This is correct, but leads to confusing
output.

=item * When calculating the probability of a sequence, there are two
serious problems.

=over

=item * There is no sequence weighting. If two sequences are almost
identical, they go into the calculation with the same weight. This is
very hard to fix without doing a clustering of sequences, but there
could be one way to cheat. Use hmmbuild, just to get the sequence
weights. Better - run a quick alignment in F<mafft> and then use our F<reduce> program
to get an even spacing over sequence space.

=item * There is no accounting for entropy in the calculation of
compatibility. This should be easy to fix. If a residue type is present 5%
of the time, it is significant when the other 95 % is conserved. If
the other 95 % of the time is scattered amongst the other amino acids,
it is not significant. Find a reference for how other methods treat
this. Explicitly, an incompatibility is interesting if the entropy is
=low. It is boring if the entropy is high.


=back

=item * When reading fasta format, walk down the start and end of the
sequences and replace gap characters with C<not_a_char>. At the moment
the number of homologues column is not meaningful with fasta format.

=item * Should we output the gaps column as it stands, or rather,
invert the sense of it. Do we want to show how many gaps there are, or
instead, how many sequences have a residue at each S<position ?>

=item * When looking for domain boundaries, there should be some way
of smoothing. N- and C- termini are not so clear.

=item * Domain boundaries. We produce a data for a plot where we mark
changes in sequence presence/absence. This will signify a domain
boundary, but it also leads to lots of peaks.

=item * DNA
We calculate entropies using log_20 because there are 20 kinds of amino acid. We should add an option like F<-4> to use log base 4 for DNA.

=item * Scaling output
If you look at the entropy within some group of proteins, the values may range from 0 to 0.2 or from 0 to 1.0. This can make it hard to put on a plot or picture. Think in the direction of quartiles. Is this point in the first or second quartile of conservation ? Now, think of a continuous version. Rank the conservation across the sites, given I<n> sites in the sequence, rank them. Each sites conservation is then given as S<I<rank / n>>.

=back

=head1 AUTHOR

Andrew Torda

=cut

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use Carp qw (cluck);
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use Seqwork qw(fasta_get_seq fasta_get_query);

# ----------------------- Constants ---------------------------------
use vars qw ($gapchar $no_char $DFLT_E_VALUE);
use vars qw ($FASTA_TYPE $BLAST_TYPE);

*gapchar      = \'-';
*no_char      = \'0';
*DFLT_E_VALUE = \0.005;
*FASTA_TYPE   = \0;  # We can force a file to be treated as
*BLAST_TYPE   = \1;  # fasta or blast format

# ----------------------- get_q_s -----------------------------------
sub get_q_s ($)
{
    my $fh = shift;
    my $line;
    my ($qstart, $qend, $qseq, $sseq);
    while ($line = <$fh>) {
        if ( ! ($line =~ m/Query/)) {
            next; }
        last;
    }
    if ( ! ($line =~ m/Query/)) {
        return undef; }
    my (@words) = split ( ' ', $line);
    $qstart = $words[1];
    $qseq =   $words[2];
    $qend =   $words[3];

    while ($line = <$fh>) {
        if ( ! ($line =~ m/^Sbjct/)) {
            next; }
        last;
    }

    if (! ($line =~ m/^Sbjct/)) {
        return undef; }
    undef (@words);
    @words = split ( ' ', $line);
    $sseq = $words[2];
    return ($qstart, $qend, $qseq, $sseq);
}

my $first = 1;
# ----------------------- get_an_align ------------------------------
sub get_an_align ($)
{
    my $fh = shift;
    my $line;
    while ( $line = <$fh>) {
        if ($line =~ m/^\>/) {
            last; } }
    if ( ! $line) {
        return (); }
    my ($score, $e_val, $frac_ident);
    while ( $line = <$fh>) {
        if ($line =~ m/Score.+Expect/) {
            $line =~ s/,/ /g;
            my @words = split (' ', $line);
            $score = $words[2];
            my $t = $words[7];   # This field is the e-value
            if ( $t =~ m/^e/ ) { # convert e-100 to 1e-100
                $t = '1' . $t; }
            $e_val = 1.0 * $t;   # Force numeric context
            last;
        }
    }
    my ($ident, $length);
    if ( ! ($line =~ m/Score.+Expect/)) {
        return (); }
    while ($line = <$fh> ) {
        if ($line =~ m/Identi.+Posi.+\%/) {
            my @words = split (' ', $line);
            ($ident, $length) = split (/\//, $words[2]);
            last;
        }
    }
    if ( ! ($line =~ /Identi.+Posi/)) {
        return (); }
    my $res_done = 0;
    my ($query, $sbjct, $start, $end);
    while ($res_done < $length) {
        my ($qstart, $qend, $qseq, $sseq) = get_q_s ($fh);
        $query .= $qseq;
        $sbjct .= $sseq;
        if ( ! $start) {
            $start = $qstart };
        $end = $qend;
        $res_done += length ($qseq);
    }
    my %align;
    $align { start }  = $start;
    $align { end }    = $end;
    $align { query }  = $query;
    $align { sbjct }  = $sbjct;
    $align { ident }  = $ident;
    $align { length } = $length;
    $align { score }  = $score;
    $align { e_val }  = $e_val;
    return ( %align );
}

# ----------------------- squash  -----------------------------------
# This is for removing query gaps from the query and subject.
# It is only useful for producing alignments which are
# entirely geared to the query.
sub squash (\$ \$)
{
    my ($q, $s) = @_;
    my $gaps = 0;
    my $l;
    if (($l = length ($$s)) != length ($$q)) {
        warn "Query and sequence strings diff lengths";
        return undef;
    }


#   The coding style: If there is a gap in query, we want to get rid
#   of it in both query and the subject sequence. This would lead to a
#   lot of shuffling (remove a character, copy, remove next, ..).
#   A faster approach is we mark each position for removal in the loop
#   below. We then use perl's fast substitution operator to get rid of
#   the magic characters.

    for (my $i = 0; $i < $l; $i++) {
        if (substr ( $$q, $i, 1) eq $gapchar) {
            substr ( $$q, $i, 1) = '9';  # Remove the position from query
            substr ( $$s, $i, 1) = '9';  # and subject
        }
    }
    $$s =~ s/9//g;
    $$q =~ s/9//g;
    return 1;
}

# ----------------------- no_char_str -------------------------------
# Return this many nothing characters
sub no_char_str ($)
{
    my $n = shift;
    my $tmp;
    for (my $i = 0; $i < $n; $i++) {
        $tmp .= $no_char; }
    return $tmp;
}

# ----------------------- merge_query -------------------------------
sub merge_query (\@ $ $ $ $ $)
{
    my ($algns, $mult_start, $mult_end, $seq, $start, $end) = @_;
#   Grow at start ?
    my ($lgrow, $rgrow);
    if ($$mult_start == -1) {
        $$algns[0] = $seq;
        $$mult_start = $start;
        $$mult_end   = $end;
        return;
    }
    $lgrow = $$mult_start - $start;

    if ($lgrow > 0) {
        $$algns[0] = no_char_str ($lgrow) . $$algns[0];
    } else {
        $lgrow = 0; }

    $rgrow = $end - $$mult_end;
    if ($rgrow > 0) {
        $$algns[0] .= no_char_str ($rgrow);
    } else {
        $rgrow = 0 ; }
#   Now our string is big enough.
    if ($lgrow && $rgrow) {
        substr ($$algns [0], 0) = $seq;
    } elsif ($lgrow) {
        substr ($$algns [0], 0, length ($seq)) = $seq;
    } elsif ($rgrow) {
        substr ($$algns [0], $start - $$mult_start) = $seq; }

    if ($start < $$mult_start) {
        $$mult_start = $start; }
    if ($end   > $$mult_end ) {
        $$mult_end   = $end; }

}

# ----------------------- add_sbjct ---------------------------------
sub add_sbjct (\@ $ $ $ $ $)
{
    my ($algns, $sbjct, $start, $end, $mult_start, $mult_len) = @_;
    my $ndx = $#$algns + 1;
    $$algns [$ndx] = no_char_str ($mult_len);
    substr ($$algns [$ndx], $start - $mult_start, length ($sbjct)) =
        $sbjct;
}

#-------------------- goto_last_round -----------------
# We look for the last occurence of a magic string which
# indicates the start of a blast iteration.
# If it is not found, we hop back to the point from which we
# began.
sub goto_last_round ($)
{
    use POSIX qw (SEEK_SET);
    my $infile = shift;
    my $start_str = 'from round [0-9]';
    my $fpos = tell ($infile);
    while (my $line = <$infile>) {
        if ( $line =~ m/$start_str/ ) {
            $fpos = tell ($infile); } }
    if (! $fpos) {
        warn "Magic marker string\n\"$start_str\"\nnot found\n";}
    # Re-position at appropriate line
    if (! seek ($infile, $fpos, SEEK_SET) ) {
        warn "seek fail on $infile";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

# ----------------------- usage   -----------------------------------
sub usage()
{
    warn "Usage\n  $0  [-a attfile [-t 5 | 6 | 20 | 21]] [-f string_in_fasta_file [-d dump_squashed_fasta_format_file ] ] [ -e min_e_val ] -m [-o outfilename] [-s offset] infilename\n";
    warn "or, better, type
    perldoc $0\n";
    return (EXIT_FAILURE);
}

# ----------------------- entropy  ----------------------------------
# We are given references to arrays for the input (alignments)
# and for arrays into which we will write the output.
# entrpy21 is the sequence entropy with a 21 symbol alphabet (including gaps)
# entrpy20  is the same with a 20 symbol alphabet (no gaps)
# entrpy6 is the 6 category version.
# gcount is the count of gaps
sub entropy ( $ $ $ $ $ $)
{
    my ($algns, $entrpy21, $entrpy20, $gcount, $compat, $seq_ndx) = @_;
    my $len = length ($$algns[0]);
    my @acount;

    my @entrpy6;

    my $j_off = ord ($gapchar);

#   Begin by counting how often each amino acid type is found at each site
#   in the alignment.
    for (my $i = 0; $i < $len; $i++) {
        $$gcount [$i] = 0;
        for ( my $seq = 0; $seq <= $#$algns; $seq++) {
            my $c = substr ($$algns[$seq], $i, 1);
            if ($c eq $no_char) {
                next; }
            my $ndx = ord ($c) - $j_off;
            $acount[$i][ $ndx ] += 1;
            if ($c eq $gapchar) {
                $$gcount[$i]++; }
        }
    }

    for ( my $i = 0; $i < $len ; $i++) {
        my $sum = 0;                          # Total number of residues at this site
        my $sum20 = 0;
        my $nc = $#{$acount[$i]};             # Number of symbols
        for ( my $j = 0; $j <= $nc; $j++) {   # Loop over symbols and count num residues
            if ($acount [$i][$j]) {           # was this symbol seen ?
                $sum += $acount [$i][$j];
                if ( $j != 0) {                      # If not a gap, add into separate
                    $sum20 += $acount [$i][$j]; }}}  # sum.


        $$entrpy21[$i] = $$entrpy20[$i] = 0.0;
        for ( my $j = 0; $j <= $nc; $j++) {
            my $i_am_not_gap = $j;
            if ( $acount [$i][$j] ) {
                my $p = $acount [$i] [$j] / $sum;
                $$entrpy21 [$i] -= $p * log ($p) ;
                if ($i_am_not_gap) {
                    my $q = $acount [$i] [$j] / $sum20;
                    $$entrpy20 [$i] -= $q * log ($q);
                }
            }
        }
        if (defined ($seq_ndx)) { # We want to calculate per site compatibility
            my $c = substr ($$algns[$seq_ndx], $i, 1); # Amino acid type at pos $i
            my $k = ord ($c) - $j_off;
            $$compat[$i] = $acount [$i] [$k] / $sum;
        }
    }
}

# ----------------------- cnvrt6   ----------------------------------
# This is a spectacularly dangerous routine. Take our multiple
# sequence alignment and reduce everything to a six letter
# alphabet.
# We cannot handle unknown residues neatly.
# "B" could be asp or asn, "Z" can be glu or gln, but these map into
# different characters in the six character alphabet.
sub cnvrt6 ($)
{
    my $algns = shift;
    my $len = length ($$algns [0]);
    for ( my $i = 0; $i < $len; $i++) {
        for ( my $seq = 0; $seq <= $#$algns; $seq++) {
            my $c = substr ($$algns[$seq], $i, 1);
            my $x = undef;
            SWITCH: {
                if ( $c eq $no_char) { $x = $c; last SWITCH ; };
                if ( $c eq $gapchar) { $x = $c; last SWITCH ; };
                if ( 'ALVIMCJ' =~ m/$c/i ) { $x = 'A'; last SWITCH}
                if ( 'FWYH'    =~ m/$c/i ) { $x = 'F'; last SWITCH}
                if ( 'STNQ'    =~ m/$c/i ) { $x = 'S'; last SWITCH}
                if ( 'KR'      =~ m/$c/i ) { $x = 'K'; last SWITCH}
                if ( 'DE'      =~ m/$c/i ) { $x = 'D'; last SWITCH}
                if ( 'GP'      =~ m/$c/i ) { $x = 'G'; last SWITCH}
                if ( $c eq 'X' || $c eq 'Z' || $c eq 'B') {
                    print STDERR "Unknown \'$c\' residue in seq number ",
                         $seq + 1, " position ", $i + 1, " In context ",
                    substr ($$algns[$seq], $i - 3, 7), "\n";
                    $x = $no_char;
                }
            }
            if ( defined ($x) ) {
                substr ($$algns [$seq], $i, 1) = $x;
                next;
            } else {
                warn "Unknown residue \'$c\'. Program broken\n";
                warn "Seq number ", $seq + 1, " position ", $i + 1, "\n";
                warn "Context ", substr ($$algns[$seq], $i - 3, 7), "\n";
                return EXIT_FAILURE;
            }
        }
    }
}


# ----------------------- do_count ----------------------------------
# We are given a file handle to write to and a two dimensional
# array with all the alignments. This routine calls entropy() to
# do the work and just writes it all out.
# If the $attfile is not null, then also write the output in a format suitable
# for a residue attribute in chimera.
# $q_seq is the query sequence. It is only known if we have fasta format input.
sub do_count ($ $ $ \@ $ $ $ $ $)
{
    my ($fh_out, $attfile, $typflag, $algns, $mult_start,
        $res_offset, $q_seq, $seq_ndx, $domain_flag) = @_;
    my $len = length ($$algns[0]);
    my @charcount;
    my @ncount;
    my $long_heading =
"# " . localtime() . "\n# Conservation calculated from blast output.
attribute: conservation
match mode: 1-to-1
recipient: residues\n";
    my $long_head_compat =
"# Compatibility of sequence with profile
attribute: compatibility
match mode: 1-to-1
recipient: residues\n";

    my $long_head_domain =
"# Domain borders\nattribute: domain_border\nmatch mode: 1-to-1\nrecipient: residues\n";
    for (my $i = 0; $i < $len; $i++) {
        for ( my $seq = 0; $seq <= $#$algns; $seq++) {
            my $c = substr ($$algns[$seq], $i, 1);
            if ($c eq $no_char) {
                next;}
            $ncount [$i]++;
        }
    }

    my (@entrpy21, @entrpy20, @gcount);
    my (@entrpy6_21, @entrpy6_20);
    my (@compat); # The compatibility of a sequence in the probabilities of neighbours

    my $cptr = undef;

    if (defined ($seq_ndx)) {     # If we have a query sequence,
        $cptr = \@compat; }       # use the compatibility array

    entropy ( $algns, \@entrpy21, \@entrpy20, \@gcount, $cptr, $seq_ndx);

    if (cnvrt6 ($algns) == EXIT_FAILURE) {
        return EXIT_FAILURE; }
    entropy ( $algns, \@entrpy6_21, \@entrpy6_20, \@gcount, undef, undef);
#   Convert to the correct log base
    {
        my $log5  = log (5.0);
        my $log6  = log (6.0);
        my $log20 = log (20.0);
        my $log21 = log (21.0);
        for ( my $i = 0; $i < $len; $i++) {
            $entrpy20[$i]   /= $log20;
            $entrpy21[$i]   /= $log21;
            $entrpy6_20[$i] /= $log5;
            $entrpy6_21[$i] /= $log6;
        }
    }
#   We are done with entropy. Now start on the likelihood of a domain boundary.
    my @domain_changes;
    if ($domain_flag) {
        do_domains ($algns, \@domain_changes); }

    print $fh_out "#res num entrpy entrpy entrpy entrpy num frac";
    if ($domain_flag) {
        print $fh_out "   dom"; }
    if (defined ($seq_ndx) ) {
        print $fh_out "    resname compat";}
    print $fh_out "\n";
    print $fh_out "#num hom     20     21   6_20   6_21 gap present\n";
#   Let us get the number of homologues ready. This is different if we are
#   working with blast or full sequences.
    my @num_hom;
    for (my $i = 0; $i < $len; $i++) {
        if ($q_seq) {  # We are working with full sequences
            $num_hom[$i] =  $#$algns + 1 - $gcount[$i];
        } else {
            $num_hom[$i] = $ncount[$i]; }
    }
    for (my $i = 0; $i < $len; $i++) {
        my $res_num = $i + $mult_start + $res_offset;
        print $fh_out $res_num,
        " $num_hom[$i] ",  sprintf ("%.5g %.5g", $entrpy20[$i], $entrpy21[$i] ),
        sprintf (" %.5g %.5g", $entrpy6_20[$i], $entrpy6_21[$i] ),
        " $gcount[$i] ",  sprintf ("%.5g", 1.0 - $gcount[$i] / $ncount[$i]) ;
        if ($domain_flag) {
            print $fh_out " ", sprintf(" %.2g", $domain_changes[$i]); }
        if (defined($seq_ndx)) {
            print  $fh_out " ", substr($q_seq, $i, 1), sprintf (" %.2g", $compat[$i]); }
        print $fh_out "\n";
    }
    if ($attfile) {
        if ( ! open (ATTFILE, ">$attfile")) {
            warn "Open failure $attfile for output: $!\n"; return EXIT_FAILURE; }

        my $type = \@entrpy20; # default column;
        if ($typflag) {
            if ($typflag == '5') {
                $type = \@entrpy6_20 }
            elsif ($typflag == '6') {
                $type = \@entrpy6_21 }
            elsif ($typflag == '20') {
                $type = \@entrpy20 }
            else {
                $type = \@entrpy21}
        }
        print ATTFILE $long_heading;

        for (my $i = 0; $i < $len; $i++) {
            my $res_num = $i + $res_offset + 1;
            printf ATTFILE ("\t:%d\t%#g\n", $res_num, $$type[$i]);
        }
        if (defined( $seq_ndx)) {
            print ATTFILE "$long_head_compat";
            for (my $i = 0; $i < $len; $i++) {
                my $res_num = $i + $res_offset + 1;
                printf ATTFILE ("\t:%d\t%#g\n", $res_num, $$cptr[$i]);
            }
        }
        if ( $domain_flag) {
            print ATTFILE "$long_head_domain";
            for (my $i = 0; $i < $len; $i++) {
                my $res_num = $i + $res_offset + 1;
                printf ATTFILE ("\t:%d\t%#g\n", $res_num, $domain_changes[$i]);
            }
        }
        close (ATTFILE);
    }
    return EXIT_SUCCESS;
}


# ----------------------- check_num ---------------------------------
# We are given a variable. If it looks like a number, return
# happy.
# If it looks like e-100, turn it into 1e-100 and return happy.
# If it does not look like a number, return undef().
sub check_num (\$)
{
    my $num = shift;
    my ($val, $unparsed) = POSIX::strtod ($$num);
    if ($unparsed != 0) {
        my $tmp = '1' . $$num;  # prepend a "1" so e-100 becomes valid
        ($val, $unparsed) = POSIX::strtod ($tmp);
        if ($unparsed != 0) {
            return undef; }
        else {
            $$num = $tmp; }
    }
    return EXIT_SUCCESS;
}

# ----------------------- fasta_get_seqline -------------------------
# We have read the title/comment. Now try to get another line of
# sequence. Return undef at the end. Originally used chomp to clean up
# the string, but this got confused after a file had been to windows and
# back (different newline terminator).
# Now use a regexp to get rid of trailing white space
sub fasta_get_seqline ($)
{
    my $fh = shift;
    my $s = '';
    my $pos = tell ($fh);
    while (my $line = <$fh>) {
        if ($line =~ /^>/) {
            seek ($fh, $pos, SEEK_SET);
            return $s;
        }
        $line =~ s/\s+//g;
        $s = "$s$line";
        $pos = tell ($fh);
    }
#   We get here either at end of file (OK) or on error.
    if (length ($s)) {
        return $s;
    } else {
        return undef;
    }
}

# ----------------------- fasta_or_blast  ---------------------------
# Look and see if the file seems to be fasta or blast format.
sub fasta_or_blast ($)
{
    my $fh = shift;
    my $pos = tell ($fh);
    my $temp = <$fh>;
    seek ($fh, $pos, SEEK_SET);
    if ($temp =~ m/^>/) {
        return $FASTA_TYPE; }
    else {
        return $BLAST_TYPE;}
}

# ----------------------- mask_squash  ------------------------------
# We are given a sequence and a mask. Copied the masked characters
# into a new sequence which we return
sub mask_squash ($ $)
{
    my ($seq_in, $mask) = @_;
    my $len = length ($seq_in);
    my $seq_out = '';
    for (my $i = 0; $i < $len; $i++) {
        if (vec ($mask, $i, 1)) {
            my $c = substr ($seq_in, $i, 1);
            $seq_out = "$seq_out$c";
        }
    }
    return $seq_out;
}

# ----------------------- add_fasta_seq -----------------------------
# We are given a sequence with its comment. Apply the mask and add it
# to the array of sequences. If $sfh is not null, it is a reference
# to a file handle where we write the squashed sequence so it can be
# read by an alignment viewer. Do this before replacing gap characters.
# There is a section in the middle where we strip leading and trailing
# gaps. I have commented it out. If we are working with full length
# sequences, this does not seem to be the correct behaviour. If working
# with partials, you might like this. Ultimately, this could be changed
# to an option, but there are already enough options.
sub add_fasta_seq (\@ \% $ $)
{
    my ($algns, $seq, $mask, $sfh) = @_;
    my $s = mask_squash ($$seq{seq}, $mask);

    if ($sfh) {
        print $sfh $$seq{cmmt}, "\n", $s, "\n"; }
    my $do_strip_leading_gaps;
    if ($do_strip_leading_gaps) {
        my $i = 0;
        while (substr($s, $i, 1) eq $gapchar) {    # Change leading gaps to
            substr($s, $i++,1) = $no_char; }       # no_char
        $i = 1;
        while (substr ($s, -$i, 1) eq $gapchar) {  # Do same for trailing
            substr ($s, -$i++, 1) = $no_char; }    # gaps (note minus sign)
    }
    my $ndx = $#$algns + 1;
    $$algns [$ndx] = $s;
}

# ----------------------- do_domains --------------------------------
# At every site in the alignment, the residue is present/absent in
# some subset of sequences. Here, we count how many changes there are
# in this subset (in or out) between adjacent positions.
# The result goes back in the @$domain_changes array. We normalise
# the result
# For debugging a bit vector, you are going to need.. p unpack ("b*", $old_vec)
# return EXIT_FAILURE/EXIT_SUCCESS
sub do_domains (\@ \@)
{
    my ($algns, $domain_changes) = @_;
    my $nseq = $#$algns + 1;
    my $len = length ($$algns[0]);
    my ($old_vec);
    vec ($old_vec, $nseq, 1) = 0;

    for (my $j = 0; $j < $nseq; $j++ ) {
        if (my $c = substr($$algns[$j], 0, 1) ne $gapchar) {
            if ($c ne $no_char) {
                vec ($old_vec, $j, 1) = 1; } } }

    $$domain_changes[$len - 1] = 0;

    for (my $i = 1; $i < $len; $i++) {
        my $next_vec;
        vec ($next_vec, $nseq, 1) = 0;
        for ( my $j = 0; $j < $nseq; $j++) {
            if (my $c = substr ($$algns[$j], $i, 1) ne $gapchar) {
                if ($c ne $no_char) {
                    vec ($next_vec, $j, 1) = 1; } } }
        my $change_vec = $old_vec ^ $next_vec; # mark changes from site to site
        my $t = unpack("%32b*", $change_vec);  # count the number of bits set
        $$domain_changes[$i-1] += $t;
        $$domain_changes[$i]   += $t;


        $old_vec = $next_vec;
    }
    my $t = ($nseq - 1) * 2;

    for ( my $j = 0; $j < $len; $j++) {
        $$domain_changes[$j] /= $t; }

    return EXIT_SUCCESS;
}

# ----------------------- get_a_seq_set_mask ------------------------
# Read up the first sequence, return it and set the mask to all bits
# on. When we are finished, put the file position back to where we started
# from.
sub get_a_seq_set_mask ($ $ $)
{
    my ($fh, $mask, $tmp_seq) = @_;
    my $pos = tell ($fh);
    if ( ! (fasta_get_seq($fh, $tmp_seq))) {
        warn ("Failed reading first seq: $!\n"); return EXIT_FAILURE; }
    if (seek ($fh, $pos, SEEK_SET) == -1) {
        warn ("Seek fail: $!\n"); return EXIT_FAILURE}
    for (my $i = 0; $i < length ($$tmp_seq{seq}); $i++) {
        vec ($$mask, $i, 1) = 1}
    return EXIT_SUCCESS;
}

# ----------------------- mymain  -----------------------------------
sub mymain ()
{
    use POSIX qw (SEEK_SET);
    use Getopt::Std;

    my %opts;
    my $max_e_val;
    my $outplotfile;
    my $attfile;
    my $errflag;
    my $fasta_name;
    my $typflag;
    my $res_offset = 0;
    my $squash_file;
    my $domain_flag;
    if ( ! getopts ('a:d:e:f:mo:s:t:', \%opts)) {
        $errflag++; }
    if ( $opts {a}) {
        $attfile = $opts {a};}
    if ( $opts {d}) {
        $squash_file = $opts {d} }
    if ( $opts {e} ) {
        $max_e_val = $opts {e}; }
    if ( $opts {f} ) {
        $fasta_name = $opts{f}; }
    if ( $opts {m} ) {
        $domain_flag = 1;}
    if ( $opts {o} ) {
        $outplotfile = $opts {o}; }
    if ( $opts {s} ) {
        $res_offset = $opts {s}; }
    if ( $opts {t} ) {
        $typflag = $opts {t}; }

    if ($typflag && ()) {
        warn "setting -t only makes sense if -a is also set. Stopping\n";
        $errflag++;
    }

    if ($squash_file) {
        if ( ! $fasta_name ) {
            warn "-d flag set, but no sequence given with -f option\n";
            $errflag++;
        }
    }

    if ( $typflag) {
        if ($typflag != '5' && $typflag != '6' && $typflag != '20' && $typflag != '21') {
            warn "-t flag must be either 5, 6, 20 or 21. Not ", $typflag, "\n";
            $errflag++;
        }
    }

    if ($errflag) {
        warn "Error parsing options.\n";
        return (usage());
    }
    undef %opts;

#   It is pretty common to type e-100, which perl does not take
#   as a number. Check for this.
    if ( defined ($max_e_val)) {
        if ( ! defined (check_num ($max_e_val))) {
            warn "Value \"$max_e_val\" for maximum e-value not number.\n";
            return (usage());
        }
    } else {
        $max_e_val = $DFLT_E_VALUE;
    }

    if ($#ARGV < 0) {
        warn "Not enough arguments\n"; return (usage()); }

    my $infile = $ARGV[0];
    if (! open (INFILE, "<$infile")) {
        print STDERR "Open fail on $infile : $!\n";
        return EXIT_FAILURE;
    }

    my ($mult_start, $mult_end) = (-1, 0);
    my @algns;

    my $q_seq;             # The query sequence that we might print out
    my $seq_ndx = undef(); # The index of the query sequence, if using one.
    my $fasta_fmt_flag = undef();
    if ( ! $fasta_name) {  # Check, it might still be in fasta format
        if ( fasta_or_blast (\*INFILE) == $FASTA_TYPE) {
            $fasta_fmt_flag = 1; } }
    if ( $fasta_name || $fasta_fmt_flag ) {
        my $sqh = undef;
        my $count;
        my $mask = '';
        my $prnt_err = 1;

        if ($fasta_name) {
            if (($seq_ndx = fasta_get_query (\*INFILE, $fasta_name, \$mask, \$q_seq, $prnt_err)) == -1) {
                close (INFILE) ; return EXIT_FAILURE; }
        } else {
            my %tmp_seq;
            if (get_a_seq_set_mask (\*INFILE, \$mask, \%tmp_seq) == EXIT_FAILURE) {
                return EXIT_FAILURE;}
            $q_seq = $tmp_seq{seq};
        }

        my $q_seq_len = length($q_seq);
        $q_seq = mask_squash ($q_seq, $mask); # get rid of gaps in the query
        my %seq;
        if ( !seek (INFILE, 0, SEEK_SET) ) {
            warn "Seek fail: $!"; return EXIT_FAILURE; }
        if ($squash_file) {
            if (! open (SQUASH, ">$squash_file")) {
                print STDERR "Open fail for output on $squash_file : $!\n";
                return EXIT_FAILURE;
            }
            $sqh = \*SQUASH;
        }
        while (fasta_get_seq (\*INFILE, \%seq)) {
            my $s_len = length ($seq{seq});
            if ($s_len != $q_seq_len) {
                print STDERR "BROKEN. Length of query sequence: $q_seq_len\n",
                "Length of this sequence: $s_len\n",
                "Working on sequence that starts\n  $seq{cmmt}\n";
                if ($squash_file) { close SQUASH};
                close INFILE;
                return EXIT_FAILURE;
            }
            add_fasta_seq(@algns, %seq, $mask, $sqh);
            $count++;
        }
        $mult_start = 1;
        print "$count sequences used\n";
    } else {                                                 # blast format
        my $count = 0;
        my $c2 = 0;
        goto_last_round (\*INFILE);
        my $pos = tell (INFILE);
        while ( my (%align) = get_an_align (\*INFILE)) {
            if (length ($align {query}) != length ($align{sbjct})) {
                warn "mismatch of query length and subject length ",
                  length ($align{query}), "  and ", length ($align{sbjct}),
                    "\n at position ", tell (INFILE), "\n";
            }
            my $query = $align {query};
            if (! squash ($query, $align {sbjct})) {
                next; }
            my ($start, $end) = ($align{start}, $align{end});
            merge_query(@algns, \$mult_start,\$mult_end, $query, $start, $end);
            $count++;
        }
        if ( !seek (INFILE, $pos, SEEK_SET) ) {
            warn "Seek fail: $!";
            close (INFILE);
            return EXIT_FAILURE;
        }
        my $mult_len = $mult_end - $mult_start + 1;
        while ( my (%align) = get_an_align (\*INFILE)) {
            if (! squash ($align{query}, $align{sbjct})) {
                next; }
            my $sbjct = $align{sbjct};
            my ($start, $end) = ($align{start}, $align{end});
            if ( $align{e_val} > $max_e_val ) {
                next; }

            add_sbjct (@algns, $sbjct, $start, $end, $mult_start, $mult_len);
            $c2++;
        }
        print "Of ", $count + 1, " alignments, ", $c2 + 1, " accepted\n";
    }
    close INFILE;  # Reading of input alignment is finished
    if ($squash_file) {
        close SQUASH; }

    my $fh = \*STDOUT;
    if ($outplotfile) {
        if (! open (OUT, ">$outplotfile")) {
            warn "Error opening $outplotfile for output plots: $!\n";
            return EXIT_FAILURE;
        }
        $fh = \*OUT;
    }

    if ($#algns == -1) {
        warn "No sequences read. Probably format is wrong (blast/fasta)\n"; return EXIT_FAILURE }
    if (do_count ($fh, $attfile, $typflag, @algns, $mult_start, $res_offset, $q_seq, $seq_ndx, $domain_flag) == EXIT_FAILURE) {
        warn "Error from do_count. Broken.\n";
        return EXIT_FAILURE;
    }
    if ($fh != \*STDOUT) {
        close $fh; }
    return EXIT_SUCCESS;
}

# ----------------------- main    -----------------------------------
exit mymain();
