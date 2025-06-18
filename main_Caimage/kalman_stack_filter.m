%kalman stack filter
function Img_p= kalman_stack_filter(Img,G, V)
    
%     Kalman filter with each frame regarded as one measurement
%     Parameters
%     ----------
%     G : scalar
%         filter gain
%     V : scalar
%         estimated variance
% 
%     Returns
%     -------
%     None.
    if nargin <2 , G = 0.8; end
    if nargin <3 , V = 0.05; end
    
    % duplicate the first N frames to avoid the initiate artifect
    N_dup = 5;
    Img1 = zeros([size(Img, 1), size(Img, 2), size(Img, 3) + N_dup]);
    Img1(:,:,1:N_dup) = Img(:,:,1:N_dup);
    Img1(:,:,N_dup+1:end) = Img;
    
    Img = Img1;
    dimz = size(Img, 3);
    dimy = size(Img, 2);
    dimx = size(Img, 1);
    Img_s = reshape(Img, [dimx*dimy, dimz]);

    %initialization
    Ik = Img_s(:, 1); % use the first image as prediction seed
    Ek = V; % use the estimated variance

    Img_p = zeros(size(Img));
    Img_p(:,:,1) = Img(:,:,1);
    % iteration
    for k = 1:dimz-1
        % correction
        Mk = Img_s(:, k+1); % current measurement
        K = Ek/(Ek+V);   %kalman gain
        Ik = G*Ik+(1-G)*Mk+K*(Mk-Ik);  % updated correction
        Ek = Ek*(1-K); % updated estimation of variance

        %prediction
        Img_p(:,:,k+1) = reshape(Ik, [dimx, dimy]);
    end
    Img_p = Img_p(:,:,N_dup + 1:end);
end