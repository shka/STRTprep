#!/usr/bin/env ruby

bam = ARGV[0]

acc2dups = Hash.new
open("| samtools view #{bam} | cut -f 1 | sort | uniq -c").each { |line|
  dups, acc = line.strip.split(/\s/)
  acc2dups[acc] = '' if dups != '1'
}

open("| bamToBed -i #{bam}").each { |line|
  cols = line.rstrip.split(/\t/)
  next if acc2dups.key?(cols[3]) || (cols[0] =~ /^RNA_SPIKE_/ && cols[5] == '+') || cols[0] =~ /^RIBO_/
  if cols[5] == '+'
    puts([cols[0], cols[1], cols[1].to_i+1, cols[3], 1, '+'].join("\t"))
  else
    puts([cols[0], cols[2].to_i-1, cols[2], cols[3], 1, '-'].join("\t"))
  end
}
