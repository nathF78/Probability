% TP1 de Probabilites : fonctions a completer et rendre sur Moodle
% Nom : Foucher
% Prénom : Nathan
% Groupe : 1SN-C

function varargout = fonctions_TP1_proba(varargin)

    switch varargin{1}     
        case 'ecriture_RVB'
            varargout{1} = feval(varargin{1},varargin{2:end});
        case {'vectorisation_par_colonne','decorrelation_colonnes'}
            [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});
        case 'calcul_parametres_correlation'
            [varargout{1},varargout{2},varargout{3}] = feval(varargin{1},varargin{2:end}); 
    end

end

% Fonction ecriture_RVB (exercice_0.m) ------------------------------------
% (Copiez le corps de la fonction ecriture_RVB du fichier du meme nom)
function image_RVB = ecriture_RVB(image_originale)
    [m,n]=size(image_originale);
    image_RVB=zeros(m/2,n/2,3);
    image_RVB(:,:,1)=image_originale(1:2:end,2:2:end);
    image_RVB(:,:,2)=(image_originale(1:2:end,1:2:end)+image_originale(2:2:end,2:2:end))/2;
    image_RVB(:,:,3)=image_originale(2:2:end,1:2:end);
end

% Fonction vectorisation_par_colonne (exercice_1.m) -----------------------
function [Vd,Vg] = vectorisation_par_colonne(I)
    Ivg=I(:,1:end-1);
    Ivd=I(:,2:end);
    Vd=Ivd(:);
    Vg=Ivg(:);
end

% Fonction calcul_parametres_correlation (exercice_1.m) -------------------
function [r,a,b] = calcul_parametres_correlation(Vd,Vg)
    n=size(Vd);
    one=ones(n(1),1);
    md=mean(Vd);
    mg=mean(Vg);
    Ed=(mean(Vd.^2-(md^2)*one))^(1/2);
    Eg=(mean(Vg.^2-(mg^2)*one))^(1/2);
    Cov=mean(Vd.*Vg-one*md*mg);
    r=Cov/Ed*Eg;
    a=Cov/Ed^2;
    b=-a*md+mg;

end

% Fonction decorrelation_colonnes (exercice_2.m) --------------------------
function [I_decorrelee,I_min] = decorrelation_colonnes(I,I_max)
    Ivg=I(:,1:end-1);
    Ivd=I(:,2:end);
    I_decorrelee=[I(:,1),(Ivd-Ivg)]; %on garde la 1er colonne
    I_min=min(min(I_decorrelee));
end



