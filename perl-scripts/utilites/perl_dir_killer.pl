#!/usr/bin/perl
use strict;
use warnings;

my $directory = '/home/psb/temporary';

opendir (DIR, $directory) or die $!;

my $i = 0;
while (my $file = readdir(DIR)) {
	if ($file ne '.' && $file ne '..') {
		print "|$file|\n";
		system "rm /home/psb/temporary/$file";
		$i = $i + 1;
		if ($i > 300) {
		#	last;
		}
	}
}

closedir(DIR);
exit(0);
