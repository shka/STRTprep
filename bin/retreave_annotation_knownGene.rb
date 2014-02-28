#!/usr/bin/env ruby

load "bin/retreave_annotation.rb"

acc2sym = Hash.new
open("| gunzip -c #{ARGV[0]} | cut -f 1,5").each { |line|
  acc, sym = line.rstrip.split(/\t/)
  acc2sym[acc] = sym.gsub(' ', '%20')
}

tbl = ARGV[1]

retreave_annotation_class1(acc2sym, tbl)
retreave_annotation_class2(acc2sym, tbl)
retreave_annotation_class3(acc2sym, tbl)
retreave_annotation_class4(acc2sym, tbl)
retreave_annotation_class5(acc2sym, tbl)
retreave_annotation_class6(acc2sym, tbl)
retreave_annotation_class7(acc2sym, tbl)
retreave_annotation_class8(acc2sym, tbl)
