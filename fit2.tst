echo 1
ve 0
gr 0
fn xygauss
gd xygauss.dat
fi 30
q
fn gauss
gd gauss.dat
fi 20
q
fn xyquad
gd xyquad.dat
li
fi
q
ve 0
fi 10
q
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
pp
gd xygauss.dat
fn xygauss
or 0 1 2 3
wt n
fi 40
q
ru dir *.dat
ad 11 4096
lf
gd lorenz.dat
fn lorenz
or 0 1 2
wt n
fi 10
q
fn poly 10
ip 1 1 1 1 1 1 1 1 1 1
md poly10.dat  0 1 0.01
gd poly10.dat
li
fi 40
q
quit
