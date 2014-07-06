function I_texture=AAM_Vector2Appearance2D(g,ObjectPixels,texturesize)

npixels=sum(ObjectPixels(:));
if(numel(g)==npixels),
    I_texture=zeros(texturesize);
    I_texture(ObjectPixels)=g;
else
    I_texturer=zeros(texturesize);
    I_textureg=zeros(texturesize);
    I_textureb=zeros(texturesize);
    I_texturer(ObjectPixels)=g(1:npixels);
    I_textureg(ObjectPixels)=g(npixels+1:2*npixels);
    I_textureb(ObjectPixels)=g(2*npixels+1:3*npixels);
    I_texture =zeros([texturesize 3]);
    I_texture(:,:,1)=I_texturer;
    I_texture(:,:,2)=I_textureg;
    I_texture(:,:,3)=I_textureb;
end

    
 