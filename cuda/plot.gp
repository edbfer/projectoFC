i = 600;
while((1)){
  set title "t=".i
  show title
  splot "dados/t".(i % 5000).".txt" using 1:2:3 with dots;
  pause 0.017;
  i = i+1;
}
