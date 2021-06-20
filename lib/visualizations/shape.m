classdef shape < handle
   properties
       vertices
       faces
       
       h % handle
   end
   
   %% Constructor
   methods (Access = public)
       function [self] = shape(vertices,faces)
           self.vertices = vertices;
           self.faces = faces;
           
           self.draw();
       end
   end
   
   %% Public Methods
   methods (Access = public)
       function [] = updateAttitude(self,rotmat)
           
           % Update the vertices:
           v_new = (rotmat*self.vertices')';
           
           % Update the handle:
           set(self.h,'Vertices',v_new);
       end
       
       function [] = draw(self)
           self.h = patch('Faces',self.faces,'Vertices',self.vertices,...
                          'FaceColor',[0.5 0.5 0.5]);
            axis equal
            grid on
            rotate3d on
       end
   end
end