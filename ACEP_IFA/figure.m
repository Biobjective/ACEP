[data] = textread('res3.txt','%f');
x = 1:1:50;
plot(x,data);
xlabel('iteration'); 
ylabel('fitness'); 
