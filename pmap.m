close;clear;

data = importdata("us70006vll_24hrRapid_14-Jan-2020.txt");
p780 = data.data(6,:);
coef = 45;

% Line
a(1) = p780(2);
b(1) = p780(1);
a(2) = p780(2)+p780(4)*coef;
b(2) = p780(1)+p780(3)*coef;
geoplot(a,b,'r')
hold

% Ellipse 
C = [a(2), b(2)];  % center 
a = p780(7)*coef*2.5;
b = p780(6)*coef*2.5;
th = linspace(0,2*pi) ; 
xe = C(1)+a*cos(th) ; 
ye = C(2)+b*sin(th) ;
geoplot(xe,ye,'r')

% Map
text(p780(2), p780(1), 'P780')
geolimits([17 19],[-68 -64])
geobasemap streets