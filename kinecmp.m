function playKine(tailX1,tailY1,tailX2,tailY2)
f = figure;
hold off;
set(f,'DoubleBuffer','on');
for j=1:10
    for i=1:size(tailX1,1)
       subplot(2,1,1), plot(tailX1(i,:),tailY2(i,:),'b');
       axis([-40 0,-20,20]);
       axis equal;
       subplot(2,1,2), plot(tailX2(i,:),tailY2(i,:),'b');
       axis([-40 0,-20,20]);
       axis equal;
       pause(.1)
    end
end