function func_displaySLM(pos, posNum, input)

all_figs = findobj(0,'type','figure');
figs2keep = all_figs(1);
delete(setdiff(all_figs, figs2keep));

FigH = figure('Name','SLMdata','IntegerHandle','on','MenuBar','none','ToolBar','none',...
    'OuterPosition',[pos(posNum,1) pos(posNum,2) pos(posNum,3) pos(posNum,4)],...
    'InnerPosition',[pos(posNum,1) pos(posNum,2) pos(posNum,3) pos(posNum,4)],...
    'WindowState','fullscreen');
imshow(input,'InitialMagnification','fit','Border','tight'); 
end




