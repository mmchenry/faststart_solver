function playKine(tailX,tailY)
f = figure;
hold off;
set(f,'DoubleBuffer','on');
for j=1:10
    for i=1:size(tailX,1)
       plot(tailX(i,:),tailY(i,:),'b');
       axis([-40 0,-20,20]);
       pause(.1)
    end
end