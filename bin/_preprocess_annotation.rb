def a2s(acc2sym, acc)
  acc2sym.key?(acc) ? acc2sym[acc] : acc
end

def preprocess_annotation_class8(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class8.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    lefts = cols[8].split(/,/)
    rights = cols[9].split(/,/)
    len = lefts.length
    if len > 1
      0.upto(len-2) { |i|
        out.puts [cols[1], rights[i], lefts[i+1], "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, cols[2]].join("\t")
      }
    end
  }
  out.close
end

def preprocess_annotation_class7(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class7.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] == cols[5] && cols[2+ofs] == cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      len = lefts.length
      if len > 1
        if cols[2] == '+'
          1.upto(len-1) { |i|
            out.puts [cols[1], lefts[i], rights[i], "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
          }
        else
          (len-1).downto(1) { |i|
            out.puts [cols[1], lefts[i], rights[i], "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
          }
        end
      end
    end
  }
  out.close
end

def preprocess_annotation_class6(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class6.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] == cols[5] && cols[2+ofs] == cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        left = lefts[0].to_i
        out.puts [cols[1], left-500 < 0 ? 0 : left-500, left, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
      else
        right = rights[-1].to_i
        out.puts [cols[1], right, right+500, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
      end
    end
  }
  out.close
end

def preprocess_annotation_class5(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class5.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] == cols[5] && cols[2+ofs] == cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        out.puts [cols[1], lefts[0], rights[0], "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
      else
        out.puts [cols[1], lefts[-1], rights[-1], "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
      end
    end
  }
  out.close
end

def preprocess_annotation_class4(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class4.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] != cols[5] && cols[2+ofs] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        cdsright = cols[6].to_i
        0.upto(lefts.length-1) { |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if left < cdsright && cdsright <= right
            out.puts [cols[1], cdsright, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
          elsif cdsright < left
            out.puts [cols[1], left, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
          end
        }
      else
        cdsleft = cols[5].to_i
        (rights.length-1).downto(0) { |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if left <= cdsleft && cdsleft < right
            out.puts [cols[1], left, cdsleft, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
          elsif right < cdsleft
            out.puts [cols[1], left, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
          end
        }
      end
    end
  }
  out.close
end

def preprocess_annotation_class3(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class3.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] != cols[5] && cols[2+ofs] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      cdsleft = cols[5].to_i
      cdsright = cols[6].to_i
      0.upto(lefts.length-1) { |i|
        left = lefts[i].to_i
        right = rights[i].to_i
        if left <= cdsleft && cdsright <= right
          out.puts [cols[1], cdsleft, cdsright, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, cols[2]].join("\t")
        elsif left <= cdsleft && cdsleft <= right
          out.puts [cols[1], cdsleft, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, cols[2]].join("\t")
        elsif left <= cdsright && cdsright <= right
          out.puts [cols[1], left, cdsright, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, cols[2]].join("\t")
        elsif cdsleft <= left && right <= cdsright
          out.puts [cols[1], left, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, cols[2]].join("\t")
        end
      }
    end
  }
  out.close
end

def preprocess_annotation_class2(acc2sym, tbl, base, ofs=1)
  fp = open("| gzip -c --best > #{base}.class2.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] != cols[5] && cols[2+ofs] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        left = lefts[0].to_i
        fp.puts [cols[1], left-500, left, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
      else
        right = rights[-1].to_i
        fp.puts [cols[1], right, right+500, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
      end
    end
  }
  fp.close
end

def preprocess_annotation_class1(acc2sym, tbl, base, ofs=1)
  out = open("| gzip -c --best > #{base}.class1.bed.gz", 'w')
  open("| gunzip -c #{tbl} | cut -f #{ofs}-").each { |line|
    cols = line.rstrip.split(/\t/)
    if cols[2+ofs] != cols[5] && cols[2+ofs] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        cdsleft = cols[5].to_i
        0.upto(lefts.length-1) { |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if right < cdsleft
              out.puts [cols[1], left, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
            elsif left < cdsleft
              out.puts [cols[1], left, cdsleft, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '+'].join("\t")
            end
          }
        else
          cdsright = cols[6].to_i
          (rights.length-1).downto(0) { |i|
            left = lefts[i].to_i
            right = rights[i].to_i
            if cdsright < left
              out.puts [cols[1], left, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
            elsif cdsright < right
              out.puts [cols[1], cdsright, right, "#{a2s(acc2sym, cols[0])}|#{cols[0]}", 0, '-'].join("\t")
            end
          }
        end
      end
    }
  out.close
end
