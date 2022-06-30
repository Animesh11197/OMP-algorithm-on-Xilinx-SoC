%% This is the script used to convert .mat datasets to create dataset files for our C code to read from.
% FILES CREATED INCLUDE:
%  data_A.dat, data_A_i.dat- conatins the real and imaginary values of Phi.
%  data_R.dat, data_R_i.dat- contains the real and imaginary values of Y.
%  K.dat - contains the sparsity for each datapoint.
%  Truth.dat - contains the golden band occupany values for all bands.

% These files contain ONLY the test set of the dataset (i.e. last 20% of the dataset)
dlmwrite("data_A.dat",real(phi(:,:,5601)),'delimiter','\n')
dlmwrite("data_R.dat",real(Y(:,:,5601)),'delimiter','\n')
dlmwrite("data_R_i.dat",imag(Y(:,:,5601)),'delimiter','\n')
dlmwrite("data_A_i.dat",imag(phi(:,:,5601)),'delimiter','\n')
dlmwrite("K.dat",sum(label_status(:,5601)),'delimiter','\n')
dlmwrite("Truth.dat",label_status(:,5601),'delimiter','\n')
for z=5602:7000
dlmwrite("data_A.dat",real(phi(:,:,z)),'-append','delimiter','\n')
dlmwrite("data_R.dat",real(Y(:,:,z)),'-append','delimiter','\n')
dlmwrite("data_R_i.dat",imag(Y(:,:,z)),'-append','delimiter','\n')
dlmwrite("data_A_i.dat",imag(phi(:,:,z)),'-append','delimiter','\n')
dlmwrite("Truth.dat",label_status(:,z),'-append','delimiter','\n')
dlmwrite("K.dat",sum(label_status(:,z)),'-append','delimiter','\n')
end