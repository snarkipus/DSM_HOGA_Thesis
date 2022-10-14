function lollipop(x,y,color,lw)
% lollipop(x,y,color,lw)	Plot lollipops (o's and sticks)
if nargin<4
    lw = 2;
	if nargin<3
	    color = '';
	end
end

h = ishold;
hold on;

%Plot circles
plot(x,y,['d' color], 'linewidth', lw);

% Make x and y row vectors, then plot as sticks
x = x(:)'; y = y(:)';
x = [x;x;nan*ones(size(x))];
y = [y;zeros(2,length(y))];
plot(x(:),y(:),color, 'linewidth', lw);

if ~h
    hold off;
end
