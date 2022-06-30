%% This function implements OMP Algorithm
% INPUT: .mat DATASETS having for different SNRs
% caution: datasets should be first loaded into the workspace before invoking this
%  script.
% Meterics of Accuracy:
% 1. Accuracy_ocupancy - predicts the probability of correctly predicting
% occupied bands
% 2. Accuracy_bandmatch - predicts the probability of correctly predicting
% each band
% 3. Accuracy - predicts the probability of correctly predicting all available bands
Yref=Y;
count=0;
predicted=0;correct_predict=0;bandmatch=0;
for z=5601:7000
    k=sum(label_status(:,z));
    Y=Yref(:,:,z);
    A=phi(:,:,z);
    [P, L] = size(A);
    [~, N] = size(Y);
    A1 = zeros(P,L);
    %% Normailzing rows of A
    for i = 1:L
        A1(:,i) = A(:,i)/norm(A(:,i));   
    end
    %% M-OMP Algorithm
    R = Y;           %% Initializing Residual
    supp = [];
    Z = zeros(L, 1);
    output=zeros(14,1);
    i = 0;
    r = [];
    while(i<k)
         for j=1:L
             Z(j) = norm( A1(:,j)'*R);
         end     
         [~, maxInd] = max(Z);
         supp = [supp maxInd];												%16-26
         As = A1(:,supp);
         x = As\R; 
         R = R - As*x;
         i = i+1;
         output(supp)=1;
    end
    for i=1:14
        if label_status(i,z)==1
            predicted=predicted+1;
            if output(i)==1
                correct_predict=correct_predict+1;
            end
        end
        if output(i)==label_status(i,z)
            bandmatch=bandmatch+1;
        end
    end
    if output==label_status(:,z)
        count=count+1;
    end
end
Accuracy_bandmatch=(bandmatch / (1400 * 14)) * 100;
Accuracy_ocupancy=(correct_predict / predicted) * 100;
Accuracy=count/1400;

