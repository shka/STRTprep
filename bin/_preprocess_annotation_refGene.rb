#!/usr/bin/env ruby

load "bin/_preprocess_annotation.rb"

acc2sym = Hash.new
open("| gunzip -c #{ARGV[0]} | gcut -f 2,13").each { |line|
  acc, sym = line.rstrip.split(/\t/)
  acc2sym[acc] = sym.gsub(' ', '%20')
}

tbl = ARGV[1]
base = ARGV[2]

preprocess_annotation_class1(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class2(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class3(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class4(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class5(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class6(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class7(acc2sym, tbl, base, ofs=2)
preprocess_annotation_class8(acc2sym, tbl, base, ofs=2)
