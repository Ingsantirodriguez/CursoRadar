function[y] = remove_outliers(x, mu, thr);
y=x; y(x<mu-thr | x>mu+thr)=[];
end