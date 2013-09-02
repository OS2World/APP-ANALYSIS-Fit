echo 1
wi 0
fn 2lorenz
gd ngauss.dat
ip 380 100 600 750 100 400 1
fi 30
q
pa -1
pr 1
pa -1
pr 2
pa -1
fn xygauss
gd xygauss.dat
fi 30
q
pa -2
pr 1
pa -1
pr 2
pa -1
fn gauss
gd gauss.dat
fi 20
q
pa -1
pr 1
pa -1
pr 2
pa -1
fn xyquad
gd xyquad.dat
li
pa -1
fi
q
pa -1
ve 0
fi 10
q
pa -1
pp
wp xyquad.a
ip 1 1 1 1 1 1
pp
rp xyquad.a
pp
ip 1 1 1 1 1 1
cp 2 5
pp
li
fi
q
pa -1
pp
gd xygauss.dat
fn xygauss
or 0 1 2 3
wt n
fi 30
q
pa -1
pl
pa -1
ru dir *.dat
ad 11 4096
lf
gd lorenz.dat
fn lorenz
or 0 1 2
wt n
fi 10

q
pa -1
fn poly 10
ip 1 1 1 1 1 1 1 1 1 1
md poly10.dat  0 1 0.01
gd poly10.dat
li
pr 1
pa -1
pr 2
pa -1
fn xyquad
gd xygauss.dat
or 0 1 2 3
wt n
wi 3 7 3 7
set zra [0:10]
li
pa -1
pr 1
pa -1
pr 2
pa -1
co 
pa -1
fi 15
q
pa -1
co
pa -1
quit
