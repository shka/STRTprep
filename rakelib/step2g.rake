##
## Step 2g - bedGraphs
##

require 'csv'

def estimate_scale(t)
  infp = open(t.source)
  tmp = infp.gets.rstrip.split("\t")
  infp.close
  return tmp[1].to_i
end

rule /\.step2g$/ => [
  ->(path){[ path.sub(/\.step2g$/, '.step2e_cnt'),
             path.sub(/^tmp\//, 'out/bam/').sub(/\.step2g$/, '.bam') ] }] do |t|
  scale = estimate_scale(t)
  sh <<EOF
bamToBed -i #{t.sources[1]} \
| gawk 'BEGIN{ FS="\t"; OFS=" " }{if($6=="+"){print $1,$2,$2+1,".",0,"+"}else{print $1,$3-1,$3,".",0,"-"}}' \
| gsort -S #{80/THREADS}% -t '\t' -k 1,1 -k 2,2n \
| guniq -c \
| gawk 'BEGIN{ FS=" "; OFS="\t" }{print $2,$3,$4,$5,$1/#{scale},$7}' \
| pigz -c > #{t.name}
EOF
end
