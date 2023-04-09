%%Intitalyze and Load Votings.

clc
clearvars

load('Pct_19Julio_Cell_SinRepetidas.mat');
load('C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\19Julio\Repetidas.mat');
load('Pct_Em_Vec_Img.mat');

%%
%Assign Image Indexes.
clear Pct_Em_Vec_Img_04Agosto

yy = 0;
Pct_Em_Vec_Img_04Agosto = [];
Pct_Em_Vec_Img_Idx = [];
for y = 1:size(Pct_Em_Vec_Img,2)

    if sum(ismember(Repetidas(:,2)-1,y)) == 0
        yy = yy + 1;
        Pct_Em_Vec_Img_04Agosto(:,yy) = Pct_Em_Vec_Img(:,y);      
        Pct_Em_Vec_Img_Idx(yy,:) = y;
        
    end
end


%%
%Compute Z-Score Filtering. 

clear MultiLabel_Cats_04Agosto  MultiLabel_Emotions  MultiLabel_Occurrence_Emotions  Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New

Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New = num2cell(zeros(size(Pct_Em_19Julio_Cell_SinRepetidas,1), size(Pct_Em_19Julio_Cell_SinRepetidas,2)));


MultiLabel_Cell = cell(1,6);
MultiLabel_Vacias = [];
MultiLabel_Cats_17Agosto = [];
MultiLabel_Size_Mat = [];

load('C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\Emotions.mat');
load('MultiLabel_Cell_13Agosto.mat');

MultiLabel_Emotions = Emotions([4,7,5,2,8,3,1,6]);

Std_Factor_Vec = [0,1/4,1/3,1/2,2/3,3/4,1,2];

B_sum = zeros(size(Pct_Em_19Julio_Cell_SinRepetidas,1),8);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    x = cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(y,1:8));
    

    [s, ind] = sort(x);
    ind = fliplr(ind);
    
    s_flip = fliplr(s);
    
    s_total = [s(1:end-1),s_flip];
    s_total = s_total(s_total>0);

   
    s_norm = normalize(s_total);
    p_twoTailed = normcdf(s_norm);
    p_twoTailed_round = round(p_twoTailed,1);

    [A_norm, B_Norm] = find(s_norm>=0);
    
    s_num = s_total(B_Norm);
    
    if max(size(s_total)) == 1

        s_num = max(x);

    end

    [A,~] = ismember(x,s_num); 
    [~,B_s] = find(A>0);

    MultiLabel_Cats = B_s;


    MultiLabel_Cats_17Agosto(y-1,1:size(MultiLabel_Cats,2)) = MultiLabel_Cats;
    MultiLabel_Size_Mat(y-1,:) = size(MultiLabel_Cats,2);
    MultiLabel_Cell{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
    MultiLabel_Cell(y,2:1+size(MultiLabel_Cats,2)) = MultiLabel_Emotions(MultiLabel_Cats);
    

    
    for yy = 1:size(MultiLabel_Cats,2)

        Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New{y,MultiLabel_Cats(yy)} = Pct_Em_19Julio_Cell_SinRepetidas{y,MultiLabel_Cats(yy)};

    end

    Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New{y,9} = sum(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New(y,1:8)));
    

end

MultiLabel_Cell{1} = {'Nombre de Imagen'};
MultiLabel_Cell(1,2:end) = {'Emociones Ganadoras'};

MultiLabel_Cell(2:end,end+1) = num2cell(MultiLabel_Size_Mat);
MultiLabel_Cell(2:end,end+1) = Pct_Em_19Julio_Cell_SinRepetidas(2:end,10);
MultiLabel_Cell(2:end,end+1) = Pct_Em_19Julio_Cell_SinRepetidas(2:end,10);



MultiLabel_Cell(1,end-2) = {'# Emociones Ganadoras'};
MultiLabel_Cell(1,end-1) = {'KAPPA Original'};
MultiLabel_Cell(1,end) = {'KAPPA New'};



%%
% Calculate and assign Fleiss' Kappa.

clear Pi K_04Agosto Pi_Sort Pi_New  A_04Agosto  KAPPA_New_Idx  Pi_Sort  K_New_2  KAPPA_New_2  KAPPA_New
clear K_New_2  K_New_3

MultiLabel_Cats_04Agosto = MultiLabel_Cats_17Agosto;
K_Sort = [4,7,5,2,8,3,1,6];

for y = 1:size(MultiLabel_Cats_04Agosto,1)

    clear K_New_2  
    Pct_Em_Vec_Img_04Agosto(:,y);
    [Pi, K_04Agosto] = KAPPA_MIO_04Agosto(Pct_Em_Vec_Img_04Agosto(:,y));
    A_04Agosto = MultiLabel_Cats_04Agosto(y,:);
    KAPPA_New_Idx = A_04Agosto(MultiLabel_Cats_04Agosto(y,:)>0);

    Pi_Sort = Pi(K_Sort);
    Pi_New(y,1:size(KAPPA_New_Idx,2)) = Pi_Sort(KAPPA_New_Idx);  
    KAPPA_New(y,:) = K_04Agosto;
    
    
    Pct_Em_Vec_Img_Sort = Pct_Em_Vec_Img_04Agosto([K_Sort,9],y);
    
    Pct_Em_Vec_Img_Sort_Mat(:,y) = Pct_Em_Vec_Img_Sort;

    K_New_2(1:size(KAPPA_New_Idx,2),1) = Pct_Em_Vec_Img_Sort(KAPPA_New_Idx);
    K_New_2(size(KAPPA_New_Idx,2)+1:size(K_Sort,2),1) = zeros(size(K_Sort,2)-size(KAPPA_New_Idx,2),1);
    
    K_New_2(:,2) = Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New{y+1,9} - K_New_2(:,1);
    [K_04Agosto_New_2] = fleiss_09Abril(K_New_2);
    KAPPA_New_2(y,:) = K_04Agosto_New_2;
    
    K_New_3(:,y) = K_New_2(:,1);

end

Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New(2:end,10) = num2cell(KAPPA_New_2);
MultiLabel_Cell(2:end,end) = num2cell(KAPPA_New_2);

%%
%Rearrange Votings after APS and Z-Score.


Pct_Em_19Julio_Cell_SinRepetidas_Original = Pct_Em_19Julio_Cell_SinRepetidas;
clear Pct_Em_19Julio_Cell_SinRepetidas
Pct_Em_19Julio_Cell_SinRepetidas = Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New;

%%
%Classify and Asign Images with Kappa Poor Agreement.
clear KAPPA_Poor_SOLOULSA A1 A11 

KAPPA_Poor_SOLOULSA = cell(1,1);
KAPPA_Grupos_KAPPA_Votantes = cell(1,1);

for y = 1:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if Pct_Em_19Julio_Cell_SinRepetidas{y,10} < 0
        KAPPA_Poor_SOLOULSA{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
        KAPPA_Poor_SOLOULSA{y,2} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
        KAPPA_Poor_SOLOULSA{y,3} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};     
    end
    
   
end


[A1,~]=find(cell2mat(KAPPA_Poor_SOLOULSA(2:end,1))<=0);
[A11,~]=find(cell2mat(KAPPA_Poor_SOLOULSA(2:end,1))==0);

KAPPA_Poor_SOLOULSA{size(Pct_Em_19Julio_Cell_SinRepetidas,1)+1,1} = size(A1,1) + size(A11,1);


KAPPA_Poor_SOLOULSA{1,1} = 'KAPPA ULSA_Cinves';
KAPPA_Poor_SOLOULSA{1,2} = 'Nombre de Imagen ULSA_Cinves';
KAPPA_Poor_SOLOULSA{1,3} = 'Participantes_ULSA_Cinves';


Sz_K_V = size(KAPPA_Grupos_KAPPA_Votantes,1);

KAPPA_Grupos_KAPPA_Votantes{1,Sz_K_V} = size(A1,1) + size(A11,1);


KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,1} = 'Grupo';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,2} = 'Mean ULSA_Cinves';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,3} = 'Std_ULSA_Cinves';


Sz_Mean_Std = size(KAPPA_Grupos_Mean_Std_KAPPA_Votantes,1);

KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,1} = 'Poor_SOLOULSA';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,2} = mean(cell2mat(KAPPA_Poor_SOLOULSA(2:end,3)));
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,3} = std(cell2mat(KAPPA_Poor_SOLOULSA(2:end,3)));


%%
%Classify and Asign Images with Kappa Slight Agreement.

clear KAPPA_Slight_SOLOULSA A1 A11
KAPPA_Slight_SOLOULSA = cell(1,1);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if (Pct_Em_19Julio_Cell_SinRepetidas{y,10} >= 0) && (Pct_Em_19Julio_Cell_SinRepetidas{y,10} <= 0.2)
        KAPPA_Slight_SOLOULSA{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
        KAPPA_Slight_SOLOULSA{y,2} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
        KAPPA_Slight_SOLOULSA{y,3} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};     
    end
    
   
end


[A1,~]=find(cell2mat(KAPPA_Slight_SOLOULSA(2:end,1))>=0);
[A11,~]=find(cell2mat(KAPPA_Slight_SOLOULSA(2:end,1))==0);

KAPPA_Slight_SOLOULSA{size(Pct_Em_19Julio_Cell_SinRepetidas,1)+1,1} = size(A1,1) + size(A11,1);



KAPPA_Slight_SOLOULSA{1,1} = 'KAPPA ULSA_Cinves';
KAPPA_Slight_SOLOULSA{1,2} = 'Nombre de Imagen ULSA_Cinves';
KAPPA_Slight_SOLOULSA{1,3} = 'Participantes_ULSA_Cinves';

Sz_K_V = size(KAPPA_Grupos_KAPPA_Votantes,1);

KAPPA_Grupos_KAPPA_Votantes{1,Sz_K_V + 1} = size(A1,1) + size(A11,1);


Sz_Mean_Std = size(KAPPA_Grupos_Mean_Std_KAPPA_Votantes,1);

KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,1} = 'Grupo';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,2} = 'Mean ULSA_Cinves';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,3} = 'Std_ULSA_Cinves';


KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,1} = 'Slight_SOLOULSA';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,2} = mean(cell2mat(KAPPA_Slight_SOLOULSA(2:end,3)));
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,3} = std(cell2mat(KAPPA_Slight_SOLOULSA(2:end,3)));


%%
%Classify and Asign Images with Kappa Fair Agreement.

clear KAPPA_Fair_SOLOULSA A1 A11 
KAPPA_Fair_SOLOULSA = cell(1,1);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if (Pct_Em_19Julio_Cell_SinRepetidas{y,10} > 0.2) && (Pct_Em_19Julio_Cell_SinRepetidas{y,10} <= 0.4)
        KAPPA_Fair_SOLOULSA{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
        KAPPA_Fair_SOLOULSA{y,2} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
        KAPPA_Fair_SOLOULSA{y,3} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};     
    end
    
   
end


[A1,~]=find(cell2mat(KAPPA_Fair_SOLOULSA(2:end,1))>=0);
[A11,~]=find(cell2mat(KAPPA_Fair_SOLOULSA(2:end,1))==0);

KAPPA_Fair_SOLOULSA{size(Pct_Em_19Julio_Cell_SinRepetidas,1)+1,1} = size(A1,1) + size(A11,1);



KAPPA_Fair_SOLOULSA{1,1} = 'KAPPA ULSA_Cinves';
KAPPA_Fair_SOLOULSA{1,2} = 'Nombre de Imagen ULSA_Cinves';
KAPPA_Fair_SOLOULSA{1,3} = 'Participantes_ULSA_Cinves';

Sz_K_V = size(KAPPA_Grupos_KAPPA_Votantes,2);

KAPPA_Grupos_KAPPA_Votantes{1,Sz_K_V + 1} = size(A1,1) + size(A11,1);



Sz_Mean_Std = size(KAPPA_Grupos_Mean_Std_KAPPA_Votantes,1);

KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,1} = 'Grupo';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,2} = 'Mean ULSA_Cinves';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,3} = 'Std_ULSA_Cinves';


KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,1} = 'Fair_SOLOULSA';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,2} = mean(cell2mat(KAPPA_Fair_SOLOULSA(2:end,3)));
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,3} = std(cell2mat(KAPPA_Fair_SOLOULSA(2:end,3)));

%%
%Classify and Asign Images with Moderate Agreement.

clear KAPPA_Moderate_SOLOULSA A1 A11 
KAPPA_Moderate_SOLOULSA = cell(1,1);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if (Pct_Em_19Julio_Cell_SinRepetidas{y,10} > 0.4) && (Pct_Em_19Julio_Cell_SinRepetidas{y,10} <= 0.6)
        KAPPA_Moderate_SOLOULSA{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
        KAPPA_Moderate_SOLOULSA{y,2} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
        KAPPA_Moderate_SOLOULSA{y,3} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};     
    end
    
   
end


[A1,~]=find(cell2mat(KAPPA_Moderate_SOLOULSA(2:end,1))>=0);
[A11,~]=find(cell2mat(KAPPA_Moderate_SOLOULSA(2:end,1))==0);

KAPPA_Moderate_SOLOULSA{size(Pct_Em_19Julio_Cell_SinRepetidas,1)+1,1} = size(A1,1) + size(A11,1);



KAPPA_Moderate_SOLOULSA{1,1} = 'KAPPA ULSA_Cinves';
KAPPA_Moderate_SOLOULSA{1,2} = 'Nombre de Imagen ULSA_Cinves';
KAPPA_Moderate_SOLOULSA{1,3} = 'Participantes_ULSA_Cinves';

Sz_K_V = size(KAPPA_Grupos_KAPPA_Votantes,2);

KAPPA_Grupos_KAPPA_Votantes{1,Sz_K_V + 1} = size(A1,1) + size(A11,1);



Sz_Mean_Std = size(KAPPA_Grupos_Mean_Std_KAPPA_Votantes,1);

KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,1} = 'Grupo';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,2} = 'Mean ULSA_Cinves';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,3} = 'Std_ULSA_Cinves';


KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,1} = 'Moderate_SOLOULSA';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,2} = mean(cell2mat(KAPPA_Moderate_SOLOULSA(2:end,3)));
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,3} = std(cell2mat(KAPPA_Moderate_SOLOULSA(2:end,3)));


%%
%Classify and Asign Images with Kappa Substantial Agreement.

clear KAPPA_Substantial_SOLOULSA A1 A11
KAPPA_Substantial_SOLOULSA = cell(1,1);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if (Pct_Em_19Julio_Cell_SinRepetidas{y,10} > 0.6) && (Pct_Em_19Julio_Cell_SinRepetidas{y,10} <= 0.8)
        KAPPA_Substantial_SOLOULSA{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
        KAPPA_Substantial_SOLOULSA{y,2} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
        KAPPA_Substantial_SOLOULSA{y,3} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};     
    end
    
   
end


[A1,~]=find(cell2mat(KAPPA_Substantial_SOLOULSA(2:end,1))>=0);
[A11,~]=find(cell2mat(KAPPA_Substantial_SOLOULSA(2:end,1))==0);

KAPPA_Substantial_SOLOULSA{size(Pct_Em_19Julio_Cell_SinRepetidas,1)+1,1} = size(A1,1) + size(A11,1);



KAPPA_Substantial_SOLOULSA{1,1} = 'KAPPA ULSA_Cinves';
KAPPA_Substantial_SOLOULSA{1,2} = 'Nombre de Imagen ULSA_Cinves';
KAPPA_Substantial_SOLOULSA{1,3} = 'Participantes_ULSA_Cinves';

Sz_K_V = size(KAPPA_Grupos_KAPPA_Votantes,2);

KAPPA_Grupos_KAPPA_Votantes{1,Sz_K_V + 1} = size(A1,1) + size(A11,1);



Sz_Mean_Std = size(KAPPA_Grupos_Mean_Std_KAPPA_Votantes,1);

KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,1} = 'Grupo';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,2} = 'Mean ULSA_Cinves';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,3} = 'Std_ULSA_Cinves';


KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,1} = 'Substantial_SOLOULSA';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,2} = mean(cell2mat(KAPPA_Substantial_SOLOULSA(2:end,3)));
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,3} = std(cell2mat(KAPPA_Substantial_SOLOULSA(2:end,3)));


%%
%Classify and Asign Images with Kappa Almost Perfect Agreement.

clear KAPPA_Almost_Perfect_SOLOULSA A1 A11
KAPPA_Almost_Perfect_SOLOULSA = cell(1,1);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if (Pct_Em_19Julio_Cell_SinRepetidas{y,10} > 0.8) && (Pct_Em_19Julio_Cell_SinRepetidas{y,10} <= 1)
        KAPPA_Almost_Perfect_SOLOULSA{y,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
        KAPPA_Almost_Perfect_SOLOULSA{y,2} = Pct_Em_19Julio_Cell_SinRepetidas{y,15};
        KAPPA_Almost_Perfect_SOLOULSA{y,3} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};     
    end
    
   
end


[A1,~]=find(cell2mat(KAPPA_Almost_Perfect_SOLOULSA(2:end,1))>=0);
[A11,~]=find(cell2mat(KAPPA_Almost_Perfect_SOLOULSA(2:end,1))==0);

KAPPA_Almost_Perfect_SOLOULSA{size(Pct_Em_19Julio_Cell_SinRepetidas,1)+1,1} = size(A1,1) + size(A11,1);


KAPPA_Almost_Perfect_SOLOULSA{1,1} = 'KAPPA ULSA_Cinves';
KAPPA_Almost_Perfect_SOLOULSA{1,2} = 'Nombre de Imagen ULSA_Cinves';
KAPPA_Almost_Perfect_SOLOULSA{1,3} = 'Participantes_ULSA_Cinves';

Sz_K_V = size(KAPPA_Grupos_KAPPA_Votantes,2);

KAPPA_Grupos_KAPPA_Votantes{1,Sz_K_V + 1} = size(A1,1) + size(A11,1);


Sz_Mean_Std = size(KAPPA_Grupos_Mean_Std_KAPPA_Votantes,1);

KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,1} = 'Grupo';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,2} = 'Mean ULSA_Cinves';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{1,3} = 'Std_ULSA_Cinves';


KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,1} = 'Almost_Perfect_SOLOULSA';
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,2} = mean(cell2mat(KAPPA_Almost_Perfect_SOLOULSA(2:end,3)));
KAPPA_Grupos_Mean_Std_KAPPA_Votantes{Sz_Mean_Std + 1,3} = std(cell2mat(KAPPA_Almost_Perfect_SOLOULSA(2:end,3)));

%%
%Calculate and Plot Number of Images per Kappa Group after 
% Z-Score Filtering and Compare Against APS's.

clear KAPPA_Groups
C={'Poor','Slight','Fair','Mod','Subs','AP'};
load('KAPPA_Groups_04Agosto.mat');

KAPPA_Groups = cell(1,1);
idx= linspace(1,size(KAPPA_Grupos_KAPPA_Votantes,2),size(KAPPA_Grupos_KAPPA_Votantes,2));


for y = 1:size(idx,2)
    A =  idx(y);
    KAPPA_Groups{y,1} = KAPPA_Grupos_KAPPA_Votantes{1,A};
   
end


figure(1)
plot(cell2mat(KAPPA_Groups_04Agosto(:,1)),'--*r')
hold on 
plot(cell2mat(KAPPA_Groups(:,1)),'-ob')
grid on
xlim([1 6])
xlabel('KAPPA Categories','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Number of Images','FontSize',14,'FontWeight','bold','Color','k')
ylim([0 700])
set(gca, 'ytick', 0:100:700,'FontSize',14);
xtickangle(0)

legend('Attentiveness Promotion','Reliability Enhancement','Location', 'Best')

set(gca, 'xtick', 1:1:6, 'xticklabel',C);



%%
%Plot Comparison of Number of Annotators per Image:
%APS vs. Z-Score Filtering Results.

figure(2)
subplot(1,2,1)
plot(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas_Original(2:end,9)),cell2mat(Pct_Em_19Julio_Cell_SinRepetidas_Original(2:end,10)),'+r')
grid on
ylim([-0.1 max(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,10)))+0.1])
set(gca, 'ytick', -0.1:0.1:max(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,10)))+0.1);
xlabel('Number of Annotators','FontSize',14,'FontWeight','bold','Color','k')

set(gca, 'ytick', 0:0.2:1.2, 'FontSize',14);
xtickangle(0)
set(gca, 'xtick', 0:20:80, 'FontSize',14);
xlim([0 85])
title('Attentiveness Promotion')


subplot(1,2,2)
plot(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,9)),cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,10)),'+b')
grid on

ylim([-0.1 max(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,10)))+0.1])
set(gca, 'ytick', -0.1:0.1:max(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,10)))+0.1);
xlabel('Number of Annotators','FontSize',14,'FontWeight','bold','Color','k')
ylabel('KAPPA','FontSize',14,'FontWeight','bold','Color','k')

set(gca, 'ytick', 0:0.2:1.2, 'FontSize',14);
set(gca, 'xtick', 0:20:80, 'FontSize',14);
xlim([0 85])
xtickangle(0)
title('Reliability Enhancement')

%%
%Compute Number of Images and Voters per Kappa Values
%for APS Results.

K_V = cell(1,1);
K_V_KAPPA = cell(1,1);
K_V_Idx = cell(1,1);

yy = 1;
for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if sum(ismember(Pct_Em_19Julio_Cell_SinRepetidas{y,9},cell2mat(K_V(:,1)))) ==0
       yy = yy + 1;
       K_V{yy,1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
    end
end

K_V{1,1} = '# Voters';

K_V(2:end,2) = sortrows(K_V(2:end,1));
K_V{1,2} = '# Voters Sorted';

for y = 2:size(K_V,1)
    
    K_V{y,3} = sum(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas(2:end,9)) == K_V{y,2});
    
    yyy = 1;
    for yy = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
        if Pct_Em_19Julio_Cell_SinRepetidas{yy,9} == cell2mat(K_V(y,2))
            yyy = yyy + 1;

            %Kappa value for each number of annotators.
            K_V_KAPPA{y,1} = K_V{y,2};
            K_V_KAPPA{y,yyy} = Pct_Em_19Julio_Cell_SinRepetidas{yy,10};
            K_V_KAPPA{1,yyy} = 'Kappa Value';

            %Image index for each K_V_Kappa value.
            K_V_Idx{y,1} = K_V{y,2};
            K_V_Idx{y,yyy} = yy;       
            K_V_Idx{1,yyy} = 'Image Index';
        end
        
    end
    K_V_KAPPA{1,1} = '# Voters';
    K_V_Idx{1,1} = '# Voters';
    
    K_V{y,4} = mean(cell2mat(K_V_KAPPA(y,2:end)));
    K_V{y,5} = nnz(cell2mat(K_V_KAPPA(y,2:end)));
    K_V{y,6} = nnz(cell2mat(K_V_Idx(y,2:end)));
    
    K_V{1,3} = 'Occurrence per # Voters';
    K_V{1,4} = 'Mean Kappa';
    K_V{1,5} = '# Images with a Kappa Value';
    K_V{1,6} = '# Images with a Kappa Value';


end
K_V_2(:,1) = sortrows(Pct_Em_19Julio_Cell_SinRepetidas(2:end,9));

%%
%Compute Number of Images and Voters per Kappa Values
%for Z-Score Filtering Results.

K_V_04Agosto = cell(1,1);
K_V_04Agosto_KAPPA = cell(1,1);
K_V_04Agosto_Idx = cell(1,1);

yy = 1;
for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas_Original,1)
    
    if sum(ismember(Pct_Em_19Julio_Cell_SinRepetidas_Original{y,9},cell2mat(K_V_04Agosto(:,1)))) ==0
       yy = yy + 1;
       K_V_04Agosto{yy,1} = Pct_Em_19Julio_Cell_SinRepetidas_Original{y,9};
    end
end

K_V_04Agosto(2:end,2) = sortrows(K_V_04Agosto(2:end,1));


for y = 2:size(K_V_04Agosto,1)
    
    K_V_04Agosto{y,3} = sum(cell2mat(Pct_Em_19Julio_Cell_SinRepetidas_Original(2:end,9)) == K_V_04Agosto{y,2});
    yyy = 1;
    for yy = 2:size(Pct_Em_19Julio_Cell_SinRepetidas_Original,1)
        if Pct_Em_19Julio_Cell_SinRepetidas_Original{yy,9} == cell2mat(K_V_04Agosto(y,2))
            yyy = yyy + 1;
            K_V_04Agosto_KAPPA{y,1} = K_V_04Agosto{y,2};
            K_V_04Agosto_KAPPA{y,yyy} = Pct_Em_19Julio_Cell_SinRepetidas_Original{yy,10};
            
            K_V_04Agosto_Idx{y,1} = K_V_04Agosto{y,2};
            K_V_04Agosto_Idx{y,yyy} = yy;           
        end
        
    end
    K_V_04Agosto{y,5} = nnz(cell2mat(K_V_04Agosto_KAPPA(y,2:end)));
    K_V_04Agosto{y,6} = nnz(cell2mat(K_V_04Agosto_Idx(y,2:end)));
    
    K_V_04Agosto{y,4} = mean(cell2mat(K_V_04Agosto_KAPPA(y,2:end)));

end

K_V_04Agosto{1,1} = '# Voters';
K_V_04Agosto{1,2} = '# Voters Sorted';
K_V_04Agosto{1,3} = 'Occurrence per # Voters';
K_V_04Agosto{1,4} = 'Mean Kappa';
K_V_04Agosto{1,5} = '# Images with a Kappa Value';
K_V_04Agosto{1,6} = '# Images with a Kappa Value';

K_V_04Agosto_2(:,1) = sortrows(Pct_Em_19Julio_Cell_SinRepetidas_Original(2:end,9));



%%
%Kappa Values for every Image within each Kappa Category
% (Reliability Evaluation).

y10 = 0;
y20 = 0;
y30 = 0;
y40 = 0;
y50 = 0;
y60 = 0;
y70 = 0;
y80 = 0;


K_V_Means_Plot = cell(1,1);
Mean_K_V_Means_Plot = cell(1,1);

K_V_Means_Plot_NumAnnot = cell(1,1);
K_V_Means_Plot_NumImg = cell(1,1);

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas,1)
    
    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} <= 10
       y10 = y10 + 1;
       K_V_Means_Plot{1,y10+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{1,y10+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{1,y10+1} = y;

    end

    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 10 && Pct_Em_19Julio_Cell_SinRepetidas{y,9} <= 20
       y20 = y20 + 1;
       K_V_Means_Plot{2,y20+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{2,y20+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{2,y20+1} = y;
    
    end
    
    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 20 && Pct_Em_19Julio_Cell_SinRepetidas{y,9} <= 30
       y30 = y30 + 1;
       K_V_Means_Plot{3,y30+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{3,y30+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{3,y30+1} = y;
   
    end
    
    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 30 && Pct_Em_19Julio_Cell_SinRepetidas{y,9} <= 40
       y40 = y40 + 1;
       K_V_Means_Plot{4,y40+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{4,y40+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{4,y40+1} = y;
    
    end

    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 40 && Pct_Em_19Julio_Cell_SinRepetidas{y,9} <= 50
       y50 = y50 + 1;
       K_V_Means_Plot{5,y50+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{5,y50+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{5,y50+1} = y;
    
    end
    
    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 50 && Pct_Em_19Julio_Cell_SinRepetidas{y,9} <= 60
       y60 = y60 + 1;
       K_V_Means_Plot{6,y60+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{6,y60+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{6,y60+1} = y;
    
    end
    
    if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 60 && Pct_Em_19Julio_Cell_SinRepetidas{y,9} <=70
       y70 = y70 + 1;
       K_V_Means_Plot{7,y70+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{7,y70+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{7,y70+1} = y;
   
    end
    
     if Pct_Em_19Julio_Cell_SinRepetidas{y,9} > 70
       y80 = y80 + 1;
       K_V_Means_Plot{8,y80+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};
       K_V_Means_Plot_NumAnnot{8,y80+1} = Pct_Em_19Julio_Cell_SinRepetidas{y,9};
       K_V_Means_Plot_NumImg{8,y80+1} = y;
   
     end   
end

for y = 1:size(K_V_Means_Plot,1)
    Mean_K_V_Means_Plot{y,2} = mean(cell2mat(K_V_Means_Plot(y,2:end)));
    
end


%%
load('C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\20Julio\Matriz_Cell_Val.mat')
load('C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\20Julio\Matriz_Cell_Test.mat')

load('Matriz_Total_Double_Originales_34320.mat')
load('Lista_Total.mat');
load('Lista_Total_Labels.mat');


%%
clear Matriz_Pct
for y = 2:size(Lista_Total_Labels,1)
   
    Matriz_Pct(y-1,:) = cell2mat(Lista_Total_Labels(y,1:8)) ./ sum(cell2mat(Lista_Total_Labels(y,1:8)));
   
    
end
%%
%Assign AP Labels per Image, Correspondingly.

clear Labels_EMOTIC_Original
vv = 0;

for y = 2:size(Lista_Total_Labels,1)
    
    vv = vv + 1;
    Labels_EMOTIC_Original{y,1} = y;
    Labels_EMOTIC_Original(y,1:11) = Lista_Total_Labels(y,1:11);
    Labels_EMOTIC_Original{y,12} = Lista_Total_Labels{y,20};


    [B,~] = max(Matriz_Pct(y-1,:));
    C = (Matriz_Pct(y-1,:) == B);
    Em = MultiLabel_Emotions(C);
    Labels_EMOTIC_Original(y,14:13+size(Em,2)) = Em;

        
end

Labels_EMOTIC_Original(1,1:11) = Lista_Total_Labels(1,1:11);
    
%%
%Sort New Labels per Image, with regard to MultiLabel_Cell Matrix.

clear Labels_EMOTIC_Ordenado

vv = 1;
for y = 2:size(MultiLabel_Cell,1)

    for yy = 2:size(Labels_EMOTIC_Original,1)
        
        if strcmp(MultiLabel_Cell{y,1},Labels_EMOTIC_Original{yy,12}) == 1
            vv = vv + 1;
            Labels_EMOTIC_Ordenado(vv,1:19) = Labels_EMOTIC_Original(yy,:);
            Labels_EMOTIC_Ordenado{vv,20} = strcmp(MultiLabel_Cell{y,1},Labels_EMOTIC_Ordenado{y,12});
        end
        
    end
    
end

Labels_EMOTIC_Ordenado(1,1:11) = Lista_Total_Labels(1,1:11);

%%
%Compute and Plot: 1) Occurrence of each Emotion
% 2) Mean Occurrence of each Emotion
% 3) Mean Kappa of Each Emotion.

clear MultiLabel_KAPPA_Emotions MultiLabel_KAPPA_Emotions_List
MultiLabel_KAPPA_Emotions = cell(9,10);
MultiLabel_KAPPA_Emotions(2:end,2:end) = num2cell(zeros(8,9));

MultiLabel_Occurrence_Emotions = num2cell(zeros(9,6));

MultiLabel_KAPPA_Emotions_List = cell(1,1);
Pct_Em_19Julio_Cell_SinRepetidas_706 = cell(1,size(Pct_Em_19Julio_Cell_SinRepetidas_Original,2));
v1 = 0;
v2 = 0;
v3 = 0;
v4 = 0;
v5 = 0;
v6 = 0;
v7 = 0;
v8 = 0;

for y = 2:size(Pct_Em_19Julio_Cell_SinRepetidas_Original,1)

        for yy = 2:9

            if Pct_Em_19Julio_Cell_SinRepetidas_Original{y,yy-1} > 0
               MultiLabel_KAPPA_Emotions{yy,3} = MultiLabel_KAPPA_Emotions{yy,3} + Pct_Em_19Julio_Cell_SinRepetidas_Original{y,10};
               MultiLabel_KAPPA_Emotions{yy,4} = MultiLabel_KAPPA_Emotions{yy,4} + 1;
               MultiLabel_Occurrence_Emotions{yy,2} = MultiLabel_Occurrence_Emotions{yy,2} + 1;
            end

            if Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New{y,yy-1} > 0
               MultiLabel_KAPPA_Emotions{yy,6} = MultiLabel_KAPPA_Emotions{yy,6} + Pct_Em_19Julio_Cell_SinRepetidas{y,10};
               MultiLabel_KAPPA_Emotions{yy,7} = MultiLabel_KAPPA_Emotions{yy,7} + 1;
               MultiLabel_Occurrence_Emotions{yy,4} = MultiLabel_Occurrence_Emotions{yy,4} + 1;
            end


            if Labels_EMOTIC_Ordenado{y,yy-1} > 0
                if isnan(Labels_EMOTIC_Ordenado{y,10}) == 0;
                    MultiLabel_KAPPA_Emotions{yy,9} = MultiLabel_KAPPA_Emotions{yy,9} + Labels_EMOTIC_Ordenado{y,10};
                    MultiLabel_KAPPA_Emotions{yy,10} = MultiLabel_KAPPA_Emotions{yy,10} + 1;
                end

               MultiLabel_Occurrence_Emotions{yy,6} = MultiLabel_Occurrence_Emotions{yy,6} + 1;
            end


            if Pct_Em_19Julio_Cell_SinRepetidas_KAPPA_New{y,yy-1} > 0

                if yy == 2
                    v1 = v1 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v1} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                end

                 if yy == 3
                    v2 = v2 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v2} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                 end     

                if yy == 4
                    v3 = v3 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v3} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                end

                 if yy == 5
                    v4 = v4 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v4} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                 end

                if yy == 6
                    v5 = v5 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v5} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                end

                 if yy == 7
                    v6 = v6 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v6} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                 end     

                if yy == 8
                    v7 = v7 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v7} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                end

                 if yy == 9
                    v8 = v8 + 1;
                   MultiLabel_KAPPA_Emotions_List{yy,v8} = Pct_Em_19Julio_Cell_SinRepetidas{y,10};

                 end  
                
            end


        end
        
    

end

MultiLabel_KAPPA_Emotions(2:end,2) = num2cell(cell2mat(MultiLabel_KAPPA_Emotions(2:end,3)) ./ cell2mat(MultiLabel_KAPPA_Emotions(2:end,4)));
MultiLabel_KAPPA_Emotions(2:end,5) = num2cell(cell2mat(MultiLabel_KAPPA_Emotions(2:end,6)) ./ cell2mat(MultiLabel_KAPPA_Emotions(2:end,7)));
MultiLabel_KAPPA_Emotions(2:end,8) = num2cell(cell2mat(MultiLabel_KAPPA_Emotions(2:end,9)) ./ cell2mat(MultiLabel_KAPPA_Emotions(2:end,10)));


MultiLabel_Occurrence_Emotions(2:end,3) = num2cell(cell2mat(MultiLabel_Occurrence_Emotions(2:end,2)) ./ (sum(cell2mat(MultiLabel_Occurrence_Emotions(2:end,2)))));
MultiLabel_Occurrence_Emotions(2:end,5) = num2cell(cell2mat(MultiLabel_Occurrence_Emotions(2:end,4)) ./ (sum(cell2mat(MultiLabel_Occurrence_Emotions(2:end,4)))));
MultiLabel_Occurrence_Emotions(2:end,7) = num2cell(cell2mat(MultiLabel_Occurrence_Emotions(2:end,6)) ./ (sum(cell2mat(MultiLabel_Occurrence_Emotions(2:end,6)))));


MultiLabel_Occurrence_Emotions(2:end,1) = MultiLabel_Emotions;
MultiLabel_Occurrence_Emotions(1,1) = {'Emotion'};
MultiLabel_Occurrence_Emotions(1,2) = {'Sum Old Occurrence'};
MultiLabel_Occurrence_Emotions(1,3) = {'Mean Old Occurrence'};

MultiLabel_Occurrence_Emotions(1,4) = {'Sum New Occurrence'};
MultiLabel_Occurrence_Emotions(1,5) = {'Mean New Occurrence'};


%%
%Plot Mean Kappa per Emotion.

figure(18)
plot(cell2mat(MultiLabel_KAPPA_Emotions(2:end,2)),'--sr','LineWidth',3)
hold on
plot(cell2mat(MultiLabel_KAPPA_Emotions(2:end,5)),'-ob','LineWidth',3)
grid on

xlabel('Emotions','FontSize',22,'FontWeight','bold','Color','k')
ylabel('Mean KAPPA','FontSize',22,'FontWeight','bold','Color','k')
legend('Attentiveness Promotion','Reliability Enhancement','Location', 'Best')
xlim([1 8])


Labels_Ingles = cell(1,8);
Labels_Ingles{1} = {'Anger/Rage'};
Labels_Ingles{2} = {'Engagement/Anticipation'};
Labels_Ingles{3} = {'Disgust/Disconnection'};
Labels_Ingles{4} = {'Fear/Worry'};
Labels_Ingles{5} = {'Joy/Affection'};
Labels_Ingles{6} = {'Sadness/Discouragement'};
Labels_Ingles{7} = {'Surprise/Amazement'};
Labels_Ingles{8} = {'Trust/Peace'};

Labels_Ingles = [Labels_Ingles{:}];
set(gca, 'xtick', 1:1:8, 'xticklabel',Labels_Ingles,'FontSize',22);
set(gca, 'ytick', 0:.1:1,'FontSize',22);

%%
%Locate and Arrange each Image's Kappa according to its 
%Number of Labels (for Reliability Enhancement).


MultiLabel_KAPPA_Plot = cell(1,1);
Mean_MLP = num2cell(zeros(1,8));
Occurrence_MLP = cell(1,1);

y1 = 0;
y2 = 0;
y3 = 0;
y4 = 0;
y5 = 0;
y6 = 0;
y7 = 0;
y8 = 0;

for y = 2:size(MultiLabel_Cell,1)
    
    
    if MultiLabel_Cell{y,end-2} == 1
        y1 = y1 + 1;
        MultiLabel_KAPPA_Plot{1,y1} = MultiLabel_Cell{y,end};
        
    end
    
    if MultiLabel_Cell{y,end-2} == 2
        y2 = y2 + 1;
        MultiLabel_KAPPA_Plot{2,y2} = MultiLabel_Cell{y,end};
        
    end
    
    if MultiLabel_Cell{y,end-2} == 3
        y3 = y3 + 1;
        MultiLabel_KAPPA_Plot{3,y3} = MultiLabel_Cell{y,end};
        
    end
    
    if MultiLabel_Cell{y,end-2} == 4
        y4 = y4 + 1;
        MultiLabel_KAPPA_Plot{4,y4} = MultiLabel_Cell{y,end};
        
    end
    if MultiLabel_Cell{y,end-2} == 5
        y5 = y5 + 1;
        MultiLabel_KAPPA_Plot{5,y5} = MultiLabel_Cell{y,end};
        
    end 
    
    if MultiLabel_Cell{y,end-2} == 6
        y6 = y6 + 1;
        MultiLabel_KAPPA_Plot{6,y6} = MultiLabel_Cell{y,end};
        
    end
    
    if MultiLabel_Cell{y,end-2} == 7
        y7 = y7 + 1;
        MultiLabel_KAPPA_Plot{7,y7} = MultiLabel_Cell{y,end};
        
    end
    
    if MultiLabel_Cell{y,end-2} == 8
        y8 = y8 + 1;
        MultiLabel_KAPPA_Plot{8,y8} = MultiLabel_Cell{y,end};
        
    end
      
end


%%
%Compute Mean Kappa for each Number of Labels per Image
%(for Reliability Enhancement).

for y = 1:size(MultiLabel_KAPPA_Plot,1)
    
    Mean_MLP{1,y} = mean(cell2mat(MultiLabel_KAPPA_Plot(y,:)));
    
end


Occurrence_MLP{1,1} = y1;
Occurrence_MLP{1,2} = y2;
Occurrence_MLP{1,3} = y3;
Occurrence_MLP{1,4} = y4;
Occurrence_MLP{1,5} = y5;
Occurrence_MLP{1,6} = y6;
Occurrence_MLP{1,7} = y7;
Occurrence_MLP{1,8} = y8;  


%%
%Locate and Arrange each Image's Kappa according to its 
%Number of Labels (for Attentiveness Promotion).

MultiLabel_KAPPA_Plot__13Agosto = cell(1,1);
Mean_MLP__13Agosto = cell(1,1);
load('MultiLabel_Cell_13Agosto.mat');

y1 = 0;
y2 = 0;
y3 = 0;
y4 = 0;
y5 = 0;
y6 = 0;
y7 = 0;
y8 = 0;

for y = 2:size(MultiLabel_Cell_13Agosto,1)
    
   
    if MultiLabel_Cell_13Agosto{y,9} == 1
        y1 = y1 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{1,y1} = MultiLabel_Cell_13Agosto{y,10};
        
    end
    
    if MultiLabel_Cell_13Agosto{y,9} == 2
        y2 = y2 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{2,y2} = MultiLabel_Cell_13Agosto{y,10};
        
    end
    
    if MultiLabel_Cell_13Agosto{y,9} == 3
        y3 = y3 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{3,y3} = MultiLabel_Cell_13Agosto{y,10};
        
    end
    
    if MultiLabel_Cell_13Agosto{y,9} == 4
        y4 = y4 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{4,y4} = MultiLabel_Cell_13Agosto{y,10};
        
    end
    if MultiLabel_Cell_13Agosto{y,9} == 5
        y5 = y5 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{5,y5} = MultiLabel_Cell_13Agosto{y,10};
        
    end 
    
    if MultiLabel_Cell_13Agosto{y,9} == 6
        y6 = y6 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{6,y6} = MultiLabel_Cell_13Agosto{y,10};
        
    end
    
    if MultiLabel_Cell_13Agosto{y,9} == 7
        y7 = y7 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{7,y7} = MultiLabel_Cell_13Agosto{y,10};
        
    end
    
    if MultiLabel_Cell_13Agosto{y,9} == 8
        y8 = y8 + 1;
        MultiLabel_KAPPA_Plot_13Agosto{8,y8} = MultiLabel_Cell_13Agosto{y,10};
        
    end
end

%%
%Compute Mean Kappa for each Number of Labels per Image
%(for Attentiveness Promotion).

for y = 1:size(MultiLabel_KAPPA_Plot_13Agosto,1)
    
    Mean_MLP_13Agosto{1,y} = mean(cell2mat(MultiLabel_KAPPA_Plot_13Agosto(y,:)));
    
end

Occurrence_MLP_13Agosto{1,1} = y1;
Occurrence_MLP_13Agosto{1,2} = y2;
Occurrence_MLP_13Agosto{1,3} = y3;
Occurrence_MLP_13Agosto{1,4} = y4;
Occurrence_MLP_13Agosto{1,5} = y5;
Occurrence_MLP_13Agosto{1,6} = y6;
Occurrence_MLP_13Agosto{1,7} = y7;
Occurrence_MLP_13Agosto{1,8} = y8;  

%%
%Plot Mean Kappa per Number of Labels (Fig. 8)
%Plot numer of Images per Number of Labels (Fig. 10).

figure(8)
plot(cell2mat(Mean_MLP_13Agosto),'--*r')
hold on
plot(cell2mat(Mean_MLP),'-ob')
grid on

xlabel('Number of Labels per Image','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Mean KAPPA ','FontSize',14,'FontWeight','bold','Color','k')
legend('Old KAPPA','New KAPPA','Location', 'Best')
xlim([1 8])
set(gca, 'xtick', 1:1:8);
ylim([0 1])
set(gca, 'ytick', 0:0.2:1);

set(gca, 'xtick', 1:1:8, 'FontSize',14);
legend('Attentiveness Promotion','Reliability Enhancement','Location', 'Best')
xtickangle(0)


figure(10)
plot(cell2mat(Occurrence_MLP_13Agosto),'--*r')
hold on
plot(cell2mat(Occurrence_MLP),'-ob')
grid on

xlabel('Number of Labels per Image','FontSize',14,'FontWeight','bold','Color','k')
ylabel('Occurrence ','FontSize',14,'FontWeight','bold','Color','k')

legend('Old KAPPA','New KAPPA','Location', 'Best')
xlim([1 8])
set(gca, 'xtick', 1:1:8, 'FontSize',14);
ylim([0 600])
set(gca, 'ytick', 0:100:600);
legend('Attentiveness Promotion','Reliability Enhancement','Location', 'Best')
xtickangle(0)

%%
%Classify each Image's Kappa according to its Kappa Group.

clear BP BP_Eval BP_Matriz_ULSA_Cinves
C={'Poor','Slight','Fair','Moderate','Substantial','Almost_Perfect'};
for y =1:size(C,2)

       BP = strcat('KAPPA_',C{y},'_SOLOULSA(1:end-1,1)');
       BP_Eval = (eval(BP));
       
       for yy = 1:size(BP_Eval,1)
           if isempty(BP_Eval{yy}) == 1
               
               BP_Matriz_ULSA_Cinves{yy,y} = 0;
               
           else
               BP_Matriz_ULSA_Cinves{yy,y} = BP_Eval{yy};
               
           end
           
       end
end

BP_Matriz_ULSA_Cinves(1,:) = C;

%%
%Spot 1 Representative Image per Kappa Category.
%(for Reliability Enhancement).

clear Idx_min Idx_med Idx_max A_min A_med_ A_max

KAPPA_Cats_Min_Med_Max = cell(7,10);

KAPPA_Cats_Min = cell(7,1);
KAPPA_Cats_Med = cell(7,1);
KAPPA_Cats_Max = cell(7,1);


Sz_Min_Med_Max = size(KAPPA_Cats_Min_Med_Max,1); 

C={'Poor','Slight','Fair','Moderate','Substantial','Almost_Perfect'};
for y =2:size(BP_Matriz_ULSA_Cinves,2)
    
    
    BP = strcat('KAPPA_',C{y},'_SOLOULSA(2:end-1,:)');
    BP_Eval = (eval(BP));
    
    BP_Ordenada = sortrows(BP_Matriz_ULSA_Cinves(2:end,y),-1);
    
    A_min = min(cell2mat(BP_Eval(2:end,1)));
    [Idx_min,C_min] = find(cell2mat(BP_Matriz_ULSA_Cinves(2:end,y)) == A_min);
    
    
    ceili= ceil(nnz(cell2mat(BP_Matriz_ULSA_Cinves(2:end,y)))/2);
    A_med = BP_Ordenada{ceili,1};
    
    [Idx_med,C_med] = find(cell2mat(BP_Matriz_ULSA_Cinves(2:end,y)) == A_med);
    
    A_max = max(cell2mat(BP_Eval(2:end,1)));
    [Idx_max,C_max] = find(cell2mat(BP_Matriz_ULSA_Cinves(2:end,y)) == A_max);
    
    Sz_Min_Med_Max = size(KAPPA_Cats_Min_Med_Max,1); 
    
    Sz_Min_Med_Max = Sz_Min_Med_Max + 1;
    
    KAPPA_Cats_Min_Med_Max(y+1,2:4) = BP_Eval(Idx_min(1),:);
    KAPPA_Cats_Min_Med_Max(y+1,5:7) = BP_Eval(Idx_med(1),:);
    KAPPA_Cats_Min_Med_Max(y+1,8:10) = BP_Eval(Idx_max(1),:);

    KAPPA_Cats_Min_Med_Max{y+1,1} = C{y};
    


    KAPPA_Cats_Min{y+1,1} = C{y};
    KAPPA_Cats_Min(y+1,2:3) = BP_Eval(Idx_min(1),1:2);
    KAPPA_Cats_Min(y+1,4:7) = MultiLabel_Cell(Idx_min(1)+1,2:5);
    KAPPA_Cats_Min{1,1} = 'Grupo';
    KAPPA_Cats_Min{1,2} = 'KAPPA_Min';
    KAPPA_Cats_Min{1,3} = 'Img_Min';
    KAPPA_Cats_Min{1,4} = 'Emotion';
    
    KAPPA_Cats_Med{y+1,1} = C{y};
    KAPPA_Cats_Med(y+1,2:3) = BP_Eval(Idx_med(1),1:2);
    KAPPA_Cats_Med(y+1,4:7) = MultiLabel_Cell(Idx_med(1)+1,2:5);
    KAPPA_Cats_Med{1,1} = 'Grupo';
    KAPPA_Cats_Med{1,2} = 'KAPPA_Med';
    KAPPA_Cats_Med{1,3} = 'Img_Med';
    KAPPA_Cats_Med{1,4} = 'Emotion';
    
    KAPPA_Cats_Max{y+1,1} = C{y};
    KAPPA_Cats_Max(y+1,2:3) = BP_Eval(Idx_max(1),1:2);
    KAPPA_Cats_Max(y+1,4:7) = MultiLabel_Cell(Idx_max(1)+1,2:5);
    KAPPA_Cats_Max{1,1} = 'Grupo';
    KAPPA_Cats_Max{1,2} = 'KAPPA_Max';
    KAPPA_Cats_Max{1,3} = 'Img_Max';
    KAPPA_Cats_Max{1,4} = 'Emotion';
end


    
    KAPPA_Cats_Min_Med_Max{1,1} = 'Grupo';
    KAPPA_Cats_Min_Med_Max{1,2} = 'KAPPA_Min';
    KAPPA_Cats_Min_Med_Max{1,3} = 'Img_Min';
    KAPPA_Cats_Min_Med_Max{1,4} = 'Numero de Participantes Min';
    
    KAPPA_Cats_Min_Med_Max{1,5} = 'KAPPA_Med';
    KAPPA_Cats_Min_Med_Max{1,6} = 'Img_Med';
    KAPPA_Cats_Min_Med_Max{1,7} = 'Numero de Participantes Med';
    
    KAPPA_Cats_Min_Med_Max{1,8} = 'KAPPA_Max';
    KAPPA_Cats_Min_Med_Max{1,9} = 'Img_Max';
    KAPPA_Cats_Min_Med_Max{1,10} = 'Numero de Participantes Max';
    
    
    KAPPA_Cats_Min_Med_Max{2,2} = [];
    KAPPA_Cats_Min_Med_Max{2,3} = [];
    KAPPA_Cats_Min_Med_Max{2,4} = [];
 
            
%%
%Load Examples of 1 Representative Image per Kappa Category
%(for Authenticity Promotion).

load('KAPPA_Cats_Min_Med_Max_Original_14Agosto.mat')


%Locate and Assign Corresponding Names of Example of Representative
%Images per Kappa Category.

KAPPA_Cats_Min_Med_Max_Vec_Original(1:5,1) = KAPPA_Cats_Min_Med_Max_Original_14Agosto(3:7,3);
KAPPA_Cats_Min_Med_Max_Vec_Original(6:10,1) = KAPPA_Cats_Min_Med_Max_Original_14Agosto(3:7,6);
KAPPA_Cats_Min_Med_Max_Vec_Original(11:16,1) = KAPPA_Cats_Min_Med_Max_Original_14Agosto(2:7,9);

Min_Med_Max_Original = cell(1,10);
Min_Med_Max_New = cell(1,9);
yyy = 0;

for y = 1:size(KAPPA_Cats_Min_Med_Max_Vec_Original,1)

    for yy = 2:size(MultiLabel_Cell,1)
    
        if strcmp(KAPPA_Cats_Min_Med_Max_Vec_Original{y,1},MultiLabel_Cell{yy,1}) == 1
            
            yyy = yyy + 1;

            Min_Med_Max_Original(yyy,:) = MultiLabel_Cell_13Agosto(yy,:);
            Min_Med_Max_New(yyy,:) = MultiLabel_Cell(yy,:);

        end

    end

end

%%
%Load Test, Train and Validation Image Names and Pixels Sums
%for our 1,120 subset.
load('nombres_Imgs_Orig_y_ULSA_Val.mat')
load('nombres_Imgs_Orig_y_ULSA_Test.mat')
load('nombres_Imgs_Orig_y_ULSA_Train.mat')

%Load Test, Train and Validation EMOTIC Original Annotations. 
load('C:\Users\carlos\Downloads\Papers_Rekognition\Bases de Datos_Tesis\2022\EMOTIC_COCO\Annotations\Annotations.mat')


%%
%Concatenate all 3 Test, Train and Valitaion Sets of our 1,120 subset
%in a Single Matrix.

%A = readtable('val_test_train.xlsx');
%tab = A;
%New_Tab = cell(1,9);

sz_val = size(nombres_Imgs_Orig_y_ULSA_Val,1);
sz_test = size(nombres_Imgs_Orig_y_ULSA_Test,1);
sz_train = size(nombres_Imgs_Orig_y_ULSA_Train,1);

nme_Imgs_Orig_ULSA_Completa(1:sz_val,:) = nombres_Imgs_Orig_y_ULSA_Val;
nme_Imgs_Orig_ULSA_Completa(sz_val+1:sz_val+sz_test,:) = nombres_Imgs_Orig_y_ULSA_Test;
nme_Imgs_Orig_ULSA_Completa(sz_val+sz_test+1:sz_val+sz_test+sz_train,:) = nombres_Imgs_Orig_y_ULSA_Train;


%Locate the 1,120 Images' Original Names
% w.r.t. the Original EMOTIC Sets.

new_nme = cell(1,4);
new_cont = 0;
for c = 1:size(nme_Imgs_Orig_ULSA_Completa,1)
    new_Sum = sum(strcmp(new_nme,nme_Imgs_Orig_ULSA_Completa{c,4}));

    if new_Sum == 0
        new_cont = new_cont + 1;
        new_nme(new_cont,:) = nme_Imgs_Orig_ULSA_Completa(c,:);

    end

end

%% 22
%Load EMOTIC Original Test, Train and Validation Annotations Merged.
load('Annotations_23554.mat');


%Create a Matrix for each of the 1,20 Images with their most Relevant
%Features w.r.t. EMOTIC Original Annotations: 
%Folder, Name, Size, ROI of Person, Labels, Gender and Age of
%Person within ROI.

New_Tab = cell(1,9);
AAA=struct2table(merge_s);
A = readtable('val_test_train.xlsx');
tab = A;

nme_Imgs_Orig_ULSA_Completa = new_nme;
for c = 1:size(nme_Imgs_Orig_ULSA_Completa,1)
    c=c
    newStr = split(nme_Imgs_Orig_ULSA_Completa{c,3},'__');
    newStrCat = strcat(newStr{1},'.jpg');
    
    newBbox2 = split(newStr{3},'.');
    newBbox = newBbox2{1};
    
    Res = 0;  
    Res2 = tab(strcmp(tab.Filename,newStrCat),:);
    Res_Anot = AAA(strcmp(AAA.filename,newStrCat),:);

    sz_D_init = 0;
    Cat_Cell = cell(1,1);
    D = 0;

    cnt = 0;

    if size(Res2,1) > 0
        for cc = 1:size(Res_Anot{1,5}{1},2)
             clear rrr
             newB = str2num(newBbox);
             Orig_Bbox = Res_Anot{1,5}{1}(cc).body_bbox;
               
            
                if isequal(newB(1),Orig_Bbox(1)) == 1 
                    Res = cell2table(cell(1,9));
                    Res{1,2}{1} = Res_Anot{1,2};
                    Res{1,3}{1} = Res_Anot{1,1};
                    Res{1,1}{1} = num2str(c);
                    %Res{1,5}{1} = newBbox;
                    Res{1,5}{1} = ['[',num2str(Res_Anot.person{1,1}(cc).body_bbox(1)),',',num2str(Res_Anot.person{1,1}(cc).body_bbox(2)),',',...
                        num2str(Res_Anot.person{1,1}(cc).body_bbox(3)),',',num2str(Res_Anot.person{1,1}(cc).body_bbox(4)),',',']'];
                    Res{1,4}{1} = ['[',num2str(Res_Anot.image_size.n_col),',',...
                                   num2str(Res_Anot.image_size.n_row),']'];

                    
                    if  isfield(Res_Anot.person{1,1}(cc),'combined_categories') == 1
                       rrr = '[';
                        for i = 1:size(Res_Anot.person{1,1}(cc).combined_categories,2)
                            if i < size(Res_Anot.person{1,1}(cc).combined_categories,2)
                                rrr=strcat(rrr,'"',Res_Anot.person{1,1}(cc).combined_categories(i),'"',',',{' '});
                            end
    
                            if i == size(Res_Anot.person{1,1}(cc).combined_categories,2)
                                rrr=strcat(rrr,'"',Res_Anot.person{1,1}(cc).combined_categories(i),'"',']');
                            end
                        end

                        Res{1,6}{1} = rrr;
                        Res{1,7}{1} = ['[',num2str(Res_Anot.person{1,1}(cc).combined_continuous.valence),',',...
                             num2str(Res_Anot.person{1,1}(cc).combined_continuous.arousal),',',...
                        num2str(Res_Anot.person{1,1}(cc).combined_continuous.dominance),']'];
                    end

                    if  isfield(Res_Anot.person{1,1}(cc),'combined_categories') == 0
                        rrr = '[';
                        for i = 1:size(Res_Anot.person{1,1}(cc).annotations_categories.categories,2)
                            if i < size(Res_Anot.person{1,1}(cc).annotations_categories.categories,2)
                                rrr=strcat(rrr,'"',Res_Anot.person{1,1}(cc).annotations_categories.categories{i},'"',',',{' '});
                            end
    
                            if i == size(Res_Anot.person{1,1}(cc).annotations_categories.categories,2)
                                rrr=strcat(rrr,'"',Res_Anot.person{1,1}(cc).annotations_categories.categories{i},'"',']');
                            end
                        end

                        Res{1,6}{1} = rrr;
                        Res{1,7}{1} = ['[',num2str(Res_Anot.person{1,1}(cc).annotations_continuous.valence),',',...
                             num2str(Res_Anot.person{1,1}(cc).annotations_continuous.arousal),',',...
                        num2str(Res_Anot.person{1,1}(cc).annotations_continuous.dominance),']'];

                    end
                    
                    Res{1,8}{1} = Res_Anot.person{1,1}(cc).gender;
                    Res{1,9}{1} = Res_Anot.person{1,1}(cc).age;
                    New_Tab(c,:) = table2cell(Res);

                end

        end

    end


    
    if size(Res2,1) == 0
        Res = cell2table(cell(1,9));
        Res(1,2:3) = Res_Anot(1,1:2);
        Res{1,2}{1} = Res_Anot{1,2};
        Res{1,3}{1} = Res_Anot{1,1};
        Res{1,1}{1} = num2str(c);
        Res{1,5}{1} = ['[',num2str(Res_Anot.person{1,1}.body_bbox(1)),',',num2str(Res_Anot.person{1,1}.body_bbox(2)),',',...
            num2str(Res_Anot.person{1,1}.body_bbox(3)),',',num2str(Res_Anot.person{1,1}.body_bbox(4)),',',']'];
        Res{1,4}{1} = ['[',num2str(Res_Anot.image_size.n_col),',',...
                       num2str(Res_Anot.image_size.n_row),']'];
        Res{1,6}{1} = ['[',Res_Anot.person{1,1}.annotations_categories.categories{1},']'];
        Res{1,7}{1} = ['[',num2str(Res_Anot.person{1,1}.annotations_continuous.valence),',',...
                       num2str(Res_Anot.person{1,1}.annotations_continuous.arousal),...
                       ',',num2str(Res_Anot.person{1,1}.annotations_continuous.dominance),']'];
        Res{1,8}{1} = Res_Anot.person{1,1}.gender;
        Res{1,9}{1} = Res_Anot.person{1,1}.age;
        New_Tab(c,:) = table2cell(Res);
    end 

end

%%
%Check for Duplicates in New_Tab Matrix.

cn = 0;
New_Tab_Sum = zeros(1,1);
for c = 1:size(New_Tab,1)
    cn = cn + 1;
    New_Tab_Sum(cn,1) = sum(strcmp(New_Tab(:,5),New_Tab{c,5}));

end
max(New_Tab_Sum)

%%
%Create a Matrix for each of the 1,120 Images with their most Relevant
%Features w.r.t. our Reliability Enhanced Annotations: 
%Folder, Name, Size, ROI of Person, Labels, Gender and Age of
%Person within ROI.

%New_Tab_ULSA = [];
New_Tab_ULSA = New_Tab;

for c = 1:size(New_Tab,1)
    c=c
    for cc = 1:size(nme_Imgs_Orig_ULSA_Completa,1)
        newStr_ULSA = split(nme_Imgs_Orig_ULSA_Completa{cc,3},'__');
        newStrCat_ULSA = strcat(newStr_ULSA{1},'.jpg');
        newStrBbox_ULSA = split(newStr_ULSA{3},'.');

        num1 = str2num(newStrBbox_ULSA{1});
        num2 = str2num(New_Tab{c,5});           
        sum(num1 == num2);

        if (strcmp(New_Tab{c,3},newStrCat_ULSA) == 1) && (sum(num1 == num2) >= 1)
            New_Tab_ULSA{c,3} = nme_Imgs_Orig_ULSA_Completa{cc,4};
            New_Tab_ULSA{c,5} = New_Tab{c,5};%newStrBbox_ULSA{1};
        end

    end

end

New_ULSA = New_Tab_ULSA;

%%
%Check for Duplicates in New_ULSA Matrix.

New_Tab_ULSA_Sum = cell(1,2);
cn = 0;
for c = 1:size(New_Tab,1)
    New_ULSA_Sum = sum(strcmp(New_Tab_ULSA(:,3),New_Tab_ULSA{c,3}));
        cn = cn + 1;
        New_Tab_ULSA_Sum{cn,1} = New_ULSA_Sum;
        New_Tab_ULSA_Sum{cn,2} = c;
        New_Tab_ULSA_Sum{cn,3} = New_Tab_ULSA{c,3};


end
 max(cell2mat(New_Tab_ULSA_Sum(:,1)))

%%
%From a Comparison between Images and Finding Corresponding Ones,
%based on ROIs:
%Create a Matrix Similiar to Previous Ones, but
%with our New Realiability Enhanced Annotations.

clear MMM

New_Tab_Labels = New_Tab;
for c = 1:size(New_Tab_ULSA,1)

          if isempty(New_Tab_ULSA{c,3}) == 0
            Res_Labels = MultiLabel_Cell(strcmp(MultiLabel_Cell(:,1),New_Tab_ULSA{c,3}),:);
            Non_Empty = find(~cellfun(@isempty,Res_Labels(1:6)));
            New_Lab = Res_Labels(Non_Empty(2:end));

            if size(New_Lab,2) == 1
                MMM = ['["',New_Lab{1}{1},'"]'];
            end

            if size(New_Lab,2) > 1
                for cc = 1:size(New_Lab,2)
    
                    if cc == 1
                        MMM = ['[','"',New_Lab{1}{1},'"'];
                    end
    
                    if (cc > 1 && cc < size(New_Lab,2))
                        MMM = [MMM,',',' ','"',New_Lab{cc}{1},'"'];
                    end
    
                     if cc == size(New_Lab,2)
                        MMM = [MMM,',','"',New_Lab{cc}{1},'"',']'];
                    end
                end

            end

            New_Tab_Labels{c,6} = MMM;

          end



end

New_Tab_Labels_Table = cell2table(New_Tab_Labels);

%writetable(New_Tab_Labels_Table, ['C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\19Julio\New_Tab_Labels_Table.csv']);


%%
%Check for Duplicates in New_Tab_Labels Matrix.

cn = 0;
New_Tab_Lables_Sum = zeros(1,1);
for c = 1:size(New_Tab,1)

    cn = cn + 1;
    New_Tab_Labels_Sum(cn,1) = sum(strcmp(New_Tab_Labels(:,3),New_Tab_Labels{1,3}));

end
max(New_Tab_Labels_Sum)

%%
%Build Matrix with 8 Vectors of 140 Images Each in order to Select,
%Further, the New RE Test, Train and Validation Sets.

rand_idx = 0;
RRR = [];
for c = 1:8

    rand_vec = (rand_idx+1:140*c);

    rand_num = randperm(140);

    rand_vec_perm = rand_vec(rand_num);
    rand_id = randperm(140);
    rand_vec_perm(rand_id);

    RRR(:,c) = rand_vec_perm(rand_id);

    rand_idx = 140*c


end

for c = 1:8
    RRR(142,c) = max(RRR(1:140,c));
    RRR(143,c) = min(RRR(1:140,c));
end

%%
%Create New Reliability Enhanced Test, Train and Validation
%Subsets' Tables. 

ix_train = 0;
ix_test = 0;
ix_val = 0;

Train_Mat = RRR(1:98,:);
Test_Mat = RRR(99:126,:);
Val_Mat = RRR(127:140,:);

Train_Vec = [];
Test_Vec = [];
Val_Vec = [];



 for c = 1:8

     Train_Vec(ix_train+1:98*c,:) = Train_Mat(:,c);
     ix_train = 98*c;

     Test_Vec(ix_test+1:28*c,:) = Test_Mat(:,c);
     ix_test = 28*c;

     Val_Vec(ix_val+1:14*c,:) = Val_Mat(:,c);
     ix_val = 14*c;

 end

 Train_perm = randperm(size(Train_Vec,1));
 Test_perm = randperm(size(Test_Vec,1));
 Val_perm = randperm(size(Val_Vec,1));

Train_Vec_perm = Train_Vec(Train_perm);
Test_Vec_perm = Test_Vec(Test_perm);
Val_Vec_perm = Val_Vec(Val_perm);

New_Tab_Labels = sortrows(New_Tab_Labels_Table,'New_Tab_Labels3','ascend');

New_Tab_Labels_Train = New_Tab_Labels(Train_Vec_perm,:);
New_Tab_Labels_Test = New_Tab_Labels(Test_Vec_perm,:);
New_Tab_Labels_Val = New_Tab_Labels(Val_Vec_perm,:);

% writetable(New_Tab_Labels_Train, ['C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\19Julio\New_Tab_Labels_Train.csv']);
% writetable(New_Tab_Labels_Test, ['C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\19Julio\New_Tab_Labels_Test.csv']);
% writetable(New_Tab_Labels_Val, ['C:\Users\carlos\Downloads\TOTAL_Matlab_10Abril\19Julio\New_Tab_Labels_Val.csv']);


%%
%Check for Images with Wrong Name.

Categories = {'Anger','Anticipation','Disgust','Fear','Joy','Sadness','Surprise','Trust'};
cuenta  = 0;
for c = 1:size(New_Tab_Labels_Val,1)
    Train_Split = split(New_Tab_Labels_Val{c,3},'_');

    if strcmp(Train_Split{1},Categories{8}) == 1
        cuenta = cuenta + 1;
    end

end

   



%%
%Check which Images Appear more than Once 
%(Same Image Containing more than One ROI).

Categories = {'Anger','Anticipation','Disgust','Fear','Joy','Sadness','Surprise','Trust'};

Cont = 0;
B_3 = New_Tab_Labels_Table;
Imgs_1120 = cell(1,9);
Idx = cell(1,5);
for c = 1:size(B_3,1)
    Sum = sum(strcmp(table2cell(New_Tab_Labels_Table(:,3)),B_3{c,3}));
       
        if Sum >= 2

              Cont = Cont + 1;
              Idx{Cont,1} = c;
              Idx{Cont,2} = Sum;
              Idx{Cont,3} = B_3{c,3};
              Idx{Cont,4} = B_3{c,5};
              Idx{Cont,5} = B_3{c,6};
                         
        end



end

%%
%Load Folder with EMOTIC Original Annotations
paths.Annotations = 'C:\Users\carlos\Downloads\Papers_Rekognition\Bases de Datos_Tesis\2022\EMOTIC_COCO\Annotations';


%Load Annotations
load(fullfile(paths.Annotations, 'Annotations.mat'))

Annotations1 = train;
Annotations2 = val;
Annotations3 = test;


%%
%Design Matrix with Train Images and Correspondign Features.

% clear New_Train
% clear New_Train_Orig
% 
% bbb = cell2table(New_Tab);
% New_Train_2 = Annotations1(:,1:size(New_Tab_Labels_Train,1));
% 
% 
% New_Train = struct('filename',table2cell(New_Tab_Labels_Train(:,3)),'folder', [],...
%     'image_size',[],'original_database',[],'person',[]);
% 
%  for c = 1:size(New_Tab_Labels_Train,1)
%     
%     bbb_Vec = bbb(strcmp(bbb.New_Tab3,New_Tab_Labels_Train{c,3}),:);
%     row_col = str2num(New_Tab_Labels_Train{c,4}{1}); 
%     row = row_col(1);
%     col = row_col(2);
%     New_Train(c).image_size.n_col =  col;
%     New_Train(c).image_size.n_row =  row;
% 
%     New_Train(c).folder = New_Tab_Labels_Train{c,2}{1};
% 
%     train_bbox = str2num(New_Tab_Labels_Train{c,5}{1});
% 
%     New_Train(c).person.body_bbox = train_bbox;
% 
% 
%     SP = split(New_Tab_Labels_Train{c,6}{1},',');
%     
%     clear train_cats
%     train_cats = cell(1,1);
%     if size(SP,1) > 1
%         for cc = 1:size(SP,1)
%             if cc < size(SP,1)
%                 train_cats{cc} = SP{cc}(3:end-1);
%             end
%     
%             if cc == size(SP,1)
%                 train_cats{cc} = SP{cc}(2:end-2);
%             end
%             
%    
%         end
%     end
% 
%     if size(SP,1) == 1
%         for cc = 1:size(SP,1)
% 
%             if cc == size(SP,1)
%                 train_cats{cc} = SP{cc}(3:end-2);
%             end
%     
%         end
%     end
% 
%     New_Train(c).person.annotations_categories.categories = train_cats;
%     
%     New_Train(c).person.annotations_contiunous.valence = [];
%     New_Train(c).person.annotations_contiunous.arousal = [];
%     New_Train(c).person.annotations_contiunous.dominance = [];
% 
%     New_Train(c).person.gender = [];
%     New_Train(c).person.age = [];
% 
% 
%     %%%%%%%%%%%%%%%%%%
% 
% 
% 
%  end

%%
%Design Matrix with Train Images and Correspondign Features.

clear New_Train_Orig2
reng = [];
aaaa = struct2table(merge_s);
idd = [];

for c = 1:size(New_Tab_Labels_Train,1)
     sz_0 = 0;
     reng = aaaa(strcmp(aaaa.filename,New_Tab_Labels_Train{c,3}),:);
     Train_Orig_box = str2num(New_Tab_Labels_Train{c,5}{1});

     for cc = 1:size(reng,1)
%          try
            idd_cnt = 0;

            clear idd
            for i = 1:size(reng.person{cc,1},2)    
             New_Train_Orig_box = reng.person{cc,1}(i).body_bbox;
    
                 if  isequal(Train_Orig_box(1:2),New_Train_Orig_box(1:2)) == 0
                     idd_cnt = idd_cnt + 1;
                     idd(idd_cnt) = i;
                 end
    
            end

             if idd_cnt > 0
                reng.person{cc,1}(idd)=[];
             end

    
%          catch 
            sz = size(reng.person{cc,1}.annotations_categories,2);
            if sz > sz_0
                idd_keep = cc;
            end
            sz_0 = sz;
    
            
     end
     New_Train_Orig2(c,:) = reng(cc,:);

end





%%
%Design Matrix with Test Images and Correspondign Features.

% clear New_Test
% New_Test_2 = Annotations2(:,1:size(New_Tab_Labels_Test,1));
% 
% New_Test = struct('filename',table2cell(New_Tab_Labels_Test(:,3)),'folder', [],...
%     'image_size',[],'original_database',[],'person',[]);
%  for c = 1:size(New_Tab_Labels_Test,1)
%     
%     row_col = str2num(New_Tab_Labels_Test{c,4}{1}); 
%     row = row_col(1);
%     col = row_col(2);
%     New_Test(c).image_size.n_col =  col;
%     New_Test(c).image_size.n_row =  row;
% 
%     test_bbox = str2num(New_Tab_Labels_Test{c,5}{1});
%     
%     New_Test(c).person.body_bbox = test_bbox;
%     New_Test(c).folder = New_Tab_Labels_Test{c,2}{1};
% 
%     clear test_cats
%     test_cats = cell(1,1);
%     SP = split(New_Tab_Labels_Test{c,6}{1},',');
%     if size(SP,1) > 1
%         for cc = 1:size(SP,1)
%             if cc < size(SP,1)
%                 test_cats{cc} = SP{cc}(3:end-1);
%             end
%     
%             if cc == size(SP,1)
%                 test_cats{cc} = SP{cc}(2:end-2);
%             end
%     
%         end
%     end
% 
%     if size(SP,1) == 1
%         for cc = 1:size(SP,1)
% 
%             if cc == size(SP,1)
%                 test_cats{cc} = SP{cc}(3:end-2);
%             end
%     
%         end
%     end
% 
%     New_Test(c).person.annotations_categories.categories = test_cats;
%     New_Test(c).person.annotations_contiunous.valence = [];
%     New_Test(c).person.annotations_contiunous.arousal = [];
%     New_Test(c).person.annotations_contiunous.dominance = [];
% 
%     New_Test(c).person.gender = [];
%     New_Test(c).person.age = [];
%  end


%%
%Design Matrix with Test Images and Correspondign Features.

clear New_Test_Orig2
reng_test = [];
aaaa = struct2table(merge_s);
idd = [];

for c = 1:size(New_Tab_Labels_Test,1)
     sz_0 = 0;
     reng_test = aaaa(strcmp(aaaa.filename,New_Tab_Labels_Test{c,3}),:);
     Test_Orig_box = str2num(New_Tab_Labels_Test{c,5}{1});

     for cc = 1:size(reng_test,1)
%          try
            idd_cnt = 0;

            clear idd
            for i = 1:size(reng_test.person{cc,1},2)
    
             New_Test_Orig_box = reng_test.person{cc,1}(i).body_bbox;
    
                 if  isequal(Test_Orig_box(1:2),New_Test_Orig_box(1:2)) == 0
                     idd_cnt = idd_cnt + 1;
                     idd(idd_cnt) = i;
                 end
    
            end
             if idd_cnt > 0
                reng_test.person{cc,1}(idd)=[];
             end

    
%          catch 
            sz = size(reng_test.person{cc,1}.annotations_categories,2);
            if sz > sz_0
                idd_keep = cc

            end
            sz_0 = sz;
    

     end
     New_Test_Orig2(c,:) = reng_test(cc,:);

end

%%
%Design Matrix with Validation Images and Correspondign Features.

% clear New_Val
% New_Val_2 = Annotations2(:,1:size(New_Tab_Labels_Val,1));
% 
% New_Val = struct('filename',table2cell(New_Tab_Labels_Val(:,3)),'folder', [],...
%     'image_size',[],'original_database',[],'person',[]);
%  for c = 1:size(New_Tab_Labels_Val,1)
%     
%     row_col = str2num(New_Tab_Labels_Val{c,4}{1}); 
%     row = row_col(1);
%     col = row_col(2);
%     New_Val(c).image_size.n_col =  col;
%     New_Val(c).image_size.n_row =  row;
% 
%     Val_bbox = str2num(New_Tab_Labels_Val{c,5}{1});
%     
%     for m = 1:size(New_Val(c).person,2)
%         New_Val(c).person(m).body_bbox =  [];
%     end
%     New_Val(c).person.body_bbox = Val_bbox;
% 
%     New_Val(c).folder = New_Tab_Labels_Val{c,2}{1};
%     
%     clear val_cats
%     val_cats = cell(1,1);
%     SP = split(New_Tab_Labels_Val{c,6}{1},',');
%     if size(SP,1) > 1
%         for cc = 1:size(SP,1)
%             if cc < size(SP,1)
%                 val_cats{cc} = SP{cc}(3:end-1);
%             end
%     
%             if cc == size(SP,1)
%                 val_cats{cc} = SP{cc}(2:end-2);
%             end
%     
%         end
%     end
% 
%     if size(SP,1) == 1
%         for cc = 1:size(SP,1)
% 
%             if cc == size(SP,1)
%                 val_cats{cc} = SP{cc}(3:end-2);
%             end
%     
%         end
%     end
% 
%     New_Val(c).person.annotations_categories.categories = val_cats;
%     
%     New_Val(c).person.annotations_contiunous.valence = [];
%     New_Val(c).person.annotations_contiunous.arousal = [];
%     New_Val(c).person.annotations_contiunous.dominance = [];
% 
%     New_Val(c).person.gender = [];
%     New_Val(c).person.age = []; 
%  end



%%
%Design Matrix with Validation Images and Correspondign Features.

clear New_Val_Orig2
reng_Val = [];
aaaa = struct2table(merge_s);
idd = [];

for c = 1:size(New_Tab_Labels_Val,1)
     sz_0 = 0;
     reng_Val = aaaa(strcmp(aaaa.filename,New_Tab_Labels_Val{c,3}),:);
     Val_Orig_box = str2num(New_Tab_Labels_Val{c,5}{1});
     for cc = 1:size(reng_Val,1)
%          try
            idd_cnt = 0;
            clear idd
            for i = 1:size(reng_Val.person{cc,1},2)
    
             New_Val_Orig_box = reng_Val.person{cc,1}(i).body_bbox;
    
                 if  isequal(Val_Orig_box(1:2),New_Val_Orig_box(1:2)) == 0
                     idd_cnt = idd_cnt + 1;
                     idd(idd_cnt) = i;
                 end
    
            end
             if idd_cnt > 0
                reng_Val.person{cc,1}(idd)=[];
             end

    
%          catch 
            sz = size(reng_Val.person{cc,1}.annotations_categories,2);
            if sz > sz_0
                idd_keep = cc;

            end
            sz_0 = sz;
    
            
     end
     New_Val_Orig2(c,:) = reng_Val(cc,:);

end


%%
%Design Matrix with Train Images and Correspondign Features.

New_Train_Orig_Mat = struct('filename',table2cell(New_Train_Orig2(:,1)),'folder', table2cell(New_Train_Orig2(:,2)),...
    'image_size',table2cell(New_Train_Orig2(:,3)),'original_database',table2cell(New_Train_Orig2(:,4)),'person',[]);
for c = 1:size(New_Train_Orig2,1)
    
    bod_box = New_Train_Orig2.person{c,1}.body_bbox;
    New_Train_Orig_Mat(c).person.body_bbox = bod_box;
    
    if  isfield(New_Train_Orig2.person{c,1},'combined_categories') == 1
        annot = New_Train_Orig2.person{c,1}.combined_categories;
        New_Train_Orig_Mat(c).person.annotations_categories.categories = annot;

        contin = New_Train_Orig2.person{c,1}.combined_continuous;
        New_Train_Orig_Mat(c).person.annotations_continuous = contin(1);
    end

    if  isfield(New_Train_Orig2.person{c,1},'combined_categories') == 0
        annot = New_Train_Orig2.person{c,1}.annotations_categories.categories;
        New_Train_Orig_Mat(c).person.annotations_categories.categories = annot;

        contin = New_Train_Orig2.person{c,1}.annotations_continuous;
        New_Train_Orig_Mat(c).person.annotations_continuous = contin(1);
    end

    gend = New_Train_Orig2.person{c,1}.gender;
    New_Train_Orig_Mat(c).person.gender = gend;

    edad = New_Train_Orig2.person{c,1}.age;
    New_Train_Orig_Mat(c).person.age = edad;

end


%%
%Design Matrix with Test Images and Correspondign Features.

New_Test_Orig_Mat = struct('filename',table2cell(New_Test_Orig2(:,1)),'folder', table2cell(New_Test_Orig2(:,2)),...
    'image_size',table2cell(New_Test_Orig2(:,3)),'original_database',table2cell(New_Test_Orig2(:,4)),'person',[]);
for c = 1:size(New_Test_Orig2,1)
    
    bod_box = New_Test_Orig2.person{c,1}.body_bbox;
    New_Test_Orig_Mat(c).person.body_bbox = bod_box;
    
    if  isfield(New_Test_Orig2.person{c,1},'combined_categories') == 1
        annot = New_Test_Orig2.person{c,1}.combined_categories;
        New_Test_Orig_Mat(c).person.annotations_categories.categories = annot;

        contin = New_Test_Orig2.person{c,1}.combined_continuous;
        New_Test_Orig_Mat(c).person.annotations_continuous = contin(1);

    end

    if  isfield(New_Test_Orig2.person{c,1},'combined_categories') == 0
        annot = New_Test_Orig2.person{c,1}.annotations_categories.categories;
        New_Test_Orig_Mat(c).person.annotations_categories.categories = annot;

        contin = New_Test_Orig2.person{c,1}.annotations_continuous;
        New_Test_Orig_Mat(c).person.annotations_continuous = contin(1);
    end

    gend = New_Test_Orig2.person{c,1}.gender;
    New_Test_Orig_Mat(c).person.gender = gend;

    edad = New_Test_Orig2.person{c,1}.age;
    New_Test_Orig_Mat(c).person.age = edad;

end

