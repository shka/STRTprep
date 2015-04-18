##
## Step 2g - bedGraphs
##

require 'csv'

def estimate_scale(t)
  libwellid2scale = Hash.new
  CSV.table(t.source).each do |row|
    libwellid2scale["#{row[:library]}.#{row[:well]}"] = row[:spikein_reads]
  end
  return libwellid2scale[t.name.pathmap('%f').pathmap('%X')]
end

rule /\.step2g$/ => [
  ->(path){ ['tmp/byGene/samples.csv',
             path.sub(/^tmp\//, 'out/bam/').sub(/\.step2g$/, '.bam')] }] do |t|
  scale = estimate_scale(t)
  sh <<EOF
bamToBed -i #{t.sources[1]} \
| gawk 'BEGIN{ FS="\t"; OFS="\t" }{if($6=="+"){print $1,$2,$2+1,".",0,"+"}else{print $1,$3-1,$3,".",0,"-"}}' \
| gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 1,1 -k 2,2n \
| guniq -c \
| gawk 'BEGIN{ FS="\t"; OFS="\t" }{print $2,$3,$4,$5,$1/#{scale},$7}' \
| pigz -c > #{t.name}
EOF
end
