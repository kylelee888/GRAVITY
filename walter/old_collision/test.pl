$n=shift;
$col=shift;
open Error, ">error-level7-$n.ax";
open Cost, ">cost-level7-$n.ax";

print Error <<EOF;
#r 0.1
#lt "n=$n"
#ly "Average relative error in force"
#lx "\\ga\\sp2\\ep"
#y l x l
#c "\\oc" cm $col
EOF

print Cost <<EOF;
#lt "n=$n"
#ly "N\\sbops\\eb / N\\sbnaive\\eb"
#lx "\\ga\\sp2\\ep"
#y l x l
#c "\\oc"
EOF

$naive=($n*($n-1));

for ($a=1; $a<25; $a*=1.1)
{
  $line=`./gravity $n 0 5 0 100 100 $a 7 | grep GRAPH`;
  print $line;
  ($dummy,$al,$err,$ex,$app)=split / /,$line;
  printf "cost ratios: %e %e\n",$ex/$naive, $app/$naive;
  print Error "$al $err\n";
  printf Cost "#cm 1\n$al %e\n#cm 2\n$al %e\n#cm 3\n$al %e\n",$ex/$naive, $app/$naive, ($ex+$app)/$naive;
}
