function I_texture=AAM_Vector2Appearance3D(g,ObjectPixels,texturesize)

npixels=sum(ObjectPixels(:)>0);
if(numel(g)==npixels),
    I_texture=zeros(texturesize);
    I_texture(ObjectPixels>0)=g;
else
    I_texturer=zeros(texturesize);
    I_textureg=zeros(texturesize);
    I_textureb=zeros(texturesize);
    I_texturer(ObjectPixels>0)=g(1:npixels);
    I_textureg(ObjectPixels>0)=g(npixels+1:2*npixels);
    I_textureb(ObjectPixels>0)=g(2*npixels+1:3*npixels);
    I_texture =zeros([texturesize 3]);
    I_texture(:,:,1)=I_texturer;
    I_texture(:,:,2)=I_textureg;
    I_texture(:,:,3)=I_textureb;
end

    
 