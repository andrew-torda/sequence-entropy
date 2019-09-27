# 22 Jan 2015
# Functions common to other files for working with sequences.
# rcsid = $Id: Seqwork.pm,v 1.3 2015/02/10 16:48:40 torda Exp $
package Seqwork;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use POSIX qw (SEEK_SET);
use strict;
use warnings;
use vars qw (@EXPORT_OK @ISA);
@ISA = qw (Exporter);
@EXPORT_OK = qw (fasta_get_seq fasta_get_query);

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

# ----------------------- fasta_get_seq  ----------------------------
# Read just one alignment from the input file.
# When there is nothing to be found, return undef;
# To flag an error, set $seq{error}
sub fasta_get_seq ($ $)
{
    my ($fh, $seq) = @_;
    my $line;
    if (! ($line = <$fh>)) {
        return undef; }    # End of file
    if ( ! ( $line =~ /^>/ )) {
        print STDERR "No fasta comment found on line \n\"$line\"\n";
        $$seq{error} = 1;
        return undef;
    }
    $line =~ s/\s+$//;
    $$seq { cmmt } = $line;  # Comment line
    if ( ! ($$seq { seq } = fasta_get_seqline($fh))) {
        print STDERR "No sequence found after ", $$seq{cmmt}, "\n";
        $$seq{error} = 1;
        return undef;
    }

    return 1;
}


# ----------------------- fasta_get_query ---------------------------
# There should be a query sequence or reference sequence in the file.
# We look for it. If it seems to occur more than once, regard it as
# an error. If we find the sequence, we are happy and use it to get
# properties such the length and a mask telling us which positions
# are actually occupied.
# While we are here, let us save the query sequence in $q_seq
# Return -1 on error and the index of the sequence if we find it.
# if $prnt_err is defined, we do print an error when the string is not
# found. Otherwise, we have been called from somebody who will take care of
# error messages.
sub fasta_get_query ($ $ $ $ $)
{
    my ($fh, $fasta_name, $mask, $q_seq, $prnt_err) = @_;
    my %seq;
    my %tmpseq;
    my $nfound = 0;
    my @found;
    my $n = 0;
    my $seq_ndx = -1;
    while (fasta_get_seq ($fh, \%tmpseq)) {
        if ($tmpseq{error}) {
            return -1; }

        if ($tmpseq{cmmt} =~ m/$fasta_name/) {
            if ($nfound == 0) {
                %seq = %tmpseq; }
            $found[$nfound++] = $tmpseq{cmmt};
            $seq_ndx = $n;
        }
        $n++;
    }
    if ($nfound == 0) {
        if ($prnt_err) {
            print STDERR "id string \"$fasta_name\" not found\n";}
        return -1;
    } elsif ($nfound > 1) {
        print STDERR "ID string not unique. Found the following occurrences\n";
        for (my $i = 0; $i < $nfound; $i++) {
            print STDERR "   $found[$i]\n"; }
        return -1;
    }

#   We have found the sequence with the relevant string and it is stored
#   in %seq.
    $$mask = '';
    my $len = length ($seq{seq});
    $$q_seq = $seq{seq};
    for (my $i = 0; $i < $len; $i++) {
        my $char = substr ($seq {seq}, $i, 1);
        if ($char eq '-') {
            vec ($$mask, $i, 1) = 0 }
        elsif ($char =~ m/[a-zA-Z]/) {
            vec ($$mask, $i, 1) = 1 }
        else {
            print STDERR "Broke char \"$char\" in string \"", $seq{seq}, "\"\n";
            return -1;
        }
    }
    return ($seq_ndx);
}
