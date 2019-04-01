% Objects and methods for simulating and visualising a Self-Propelled Rod
% (SPR) model based on electrostatic-like interactions between rods
% composed of point-like segments.
% 
% For full details, see:
%
% Wensink, H. H., & Löwen, H. (2012). Emergent states in dense systems of 
% active rods: From swarming to turbulence. Journal of Physics Condensed 
% Matter, 24(46). https://doi.org/10.1088/0953-8984/24/46/464130
%
% Author: Oliver J. Meacock

classdef WensinkField
    properties
        xWidth %Width of field
        yHeight %Height of field
        
        xCells %x positions of centers of mass of cells in field
        yCells %Likewise for y positions
        thetCells %Likewise for orientations of cells in xy-plane
        aCells %Likewise for aspect ratios of cells
        nCells %Likewise for number of segments in cells
        uCells %Likewise for orientation vector of cells
        lCells %Likewise for separation between Yukawa segments of cells.
        fCells %Likewise for the extra force added in the direction of motion (The 'self propelled' bit)
        cCells %Likewise for color vector of cells.
        
        cellDists %Distance between all cells. Indexed in same way as obj.cells.
        distThresh %Distance between cells, beyond which interactions are too weak to be relevant. 
        U0 %Potential Amplitude
        lam %Screening length of Yukawa segements
    end
    methods
        function obj = WensinkField(xWidth,yHeight,U0,lam)
            if isnumeric(xWidth) && xWidth > 0 && isnumeric(yHeight) && yHeight > 0
                obj.xWidth = xWidth;
                obj.yHeight= yHeight;
            else
                error('Input arguments to WensinkField are not valid');
            end
            
            obj.U0 = U0;
            obj.lam = lam;
        end
        
        function obj = populateField(obj,cellSettingsType,cellSettings,areaFrac)
            %Populates the field with a number of randomly positioned cells of the given force generation and a number of static 'barrier' rods.   
            
            switch cellSettingsType
                case 'Mapped' %i.e. mapped from orientation field
                    if isinf(areaFrac) %Indicates a trial field to establish the amount of area that will be covered by rods
                        obj.xCells = 0.5*ones(cellSettings.noCells,1);
                        obj.yCells = 0.5*ones(cellSettings.noCells,1);
                        obj.aCells = cellSettings.a * ones(cellSettings.noCells,1);
                        obj.thetCells = zeros(cellSettings.noCells,1);
                    else
                        obj = obj.mapOrientations(cellSettings);
                        obj.fCells = cellSettings.f * ones(cellSettings.noCells,1);
                    end
            end
            
            [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
            obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
        end
        
        function obj = mapOrientations(obj,cS)
            for i = 1:cS.noCells
                obj.xCells(i,1) = (rand(1) * (obj.xWidth - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                obj.yCells(i,1) = (rand(1) * (obj.yHeight - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                obj.thetCells(i,1) = ((rand(1)<0.5)*pi) - cS.oriField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling));
                obj.aCells(i,1) = cS.a;
                obj.cCells(i,:) = transpose([1;1;1] - (double(squeeze(cS.colourField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling),:)))/255));
                
                %Calculate derived values
                [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
                obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
                
                obj = obj.calcDistThresh();
                obj = obj.calcDistMat();
                
                [crossCell1,~] = obj.findCrossingCells(1:length(obj.aCells),1:length(obj.aCells));
                
                while ~isempty(crossCell1)
                    obj.xCells(i,1) = (rand(1) * (obj.xWidth - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                    obj.yCells(i,1) = (rand(1) * (obj.yHeight - (2*ceil(cS.fieldScaling)+1))) + 1 + ceil(cS.fieldScaling);
                    obj.thetCells(i,1) = ((rand(1)<0.5)*pi) - cS.oriField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling));
                    obj.aCells(i,1) = cS.a;
                    obj.cCells(i,:) = transpose([1;1;1] - (double(squeeze(cS.colourField(round(obj.yCells(i) / cS.fieldScaling),round(obj.xCells(i) / cS.fieldScaling),:)))/255));
                    
                    [obj.nCells,obj.lCells] = calculateSegmentNumberLength(obj.aCells,obj.lam);
                    obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
                    
                    obj = obj.calcDistMat();
                    
                    [crossCell1,~] = obj.findCrossingCells(1:length(obj.aCells),1:length(obj.aCells));
                end
                
                progressbar(i/cS.noCells);
            end
            progressbar(1);
        end
        
        function [crossCell1,crossCell2] = findCrossingCells(obj,candList1,candList2)
            
            closeMat = obj.cellDists < max(obj.lCells.*(obj.nCells+1)); %All cells whose centroids are close enough that they could be crossing
            closeMat(logical(eye(size(closeMat)))) = 0; %Can't cross self.
            closeMat = closeMat(candList1,candList2);
            
            xs = obj.xCells;
            ys = obj.yCells;
            ls = obj.lCells;
            ns = obj.nCells;
            thets = obj.thetCells;
            
            %Prepare lists of cell start and end points
            starts1 = [xs(candList1) + ls(candList1).*(ns(candList1)-1).*cos(thets(candList1))/2,ys(candList1) + ls(candList1).*(ns(candList1)-1).*sin(thets(candList1))/2];
            ends1 = [xs(candList1) - ls(candList1).*(ns(candList1)-1).*cos(thets(candList1))/2,ys(candList1) - ls(candList1).*(ns(candList1)-1).*sin(thets(candList1))/2];
            
            starts2 = [xs(candList2) + ls(candList2).*(ns(candList2)-1).*cos(thets(candList2))/2,ys(candList2) + ls(candList2).*(ns(candList2)-1).*sin(thets(candList2))/2];
            ends2 = [xs(candList2) - ls(candList2).*(ns(candList2)-1).*cos(thets(candList2))/2,ys(candList2) - ls(candList2).*(ns(candList2)-1).*sin(thets(candList2))/2];
            
            lines1 = [starts1,ends1];
            lines2 = [starts2,ends2];
            
            crossMatStruct = lineSegmentIntersectPeriodic(lines1,lines2,closeMat,obj.xWidth/2,obj.yHeight/2,obj.xWidth,obj.yHeight);
            
            crossMat = crossMatStruct.intAdjacencyMatrix;
            
            [crossCell1,crossCell2] = ind2sub(size(crossMat),find(crossMat));
        end
        
        function areaFrac = getAreaFraction(obj) %Calculates the area fraction occupied by rods.            
            %Calculate rods as sperocylinders
            arRods = (obj.lam^2 * (obj.aCells - 1)) + (pi * (obj.lam/2)^2);
            totArea = obj.xWidth*obj.yHeight;
            areaFrac = sum(arRods)/totArea;
        end
        
        function obj = calcDistMat(obj)
            obj.cellDists = calcGriddedDistMat(obj);
        end
        
        function obj = calcDistThresh(obj)
            %Calculates the distance threshold below which interactions will be ignored. Based on the longest cell.
            foo = ((obj.nCells - 1) .* obj.lCells);
            maxLen = max(foo); %Length of the longest cell
            
            obj.distThresh = maxLen + obj.lam + log(obj.U0); %Use of log(U0) here is somewhat justified by the exponential drop off in repelling strength. But not terribly. Use cautiously.
        end
        
        function [drdt,dthetadt] = calcVelocities(obj,f0,compiled)
            %Calculates the rate of translation and rotation for all cells
            includeMat = obj.cellDists < obj.distThresh;
            includeMat(logical(diag(ones(length(includeMat),1)))) = 0;
            includeMat = includeMat(1:length(obj.nCells),:); %Alpha rods should include all cells, beta rods all cells and all barrier rods
             
            [fT,fR] = calcFrictionTensors(obj.aCells,obj.uCells,f0);
            [fPar,~,~] = calcGeomFactors(obj.aCells);
            
            drdt = zeros(size(obj.xCells,1),2);
            dthetadt = zeros(size(obj.xCells));
            
            boundX = obj.xWidth/2;
            boundY = obj.yHeight/2;
            Height = obj.yHeight;
            Width = obj.xWidth;
            lam = obj.lam;
            U0 = obj.U0;
            
            for i = 1:length(obj.nCells) %The index of the cell alpha.
                %Get indices of cells that this cell (alpha) interacts with
                betInds = find(includeMat(i,:));
                xs = obj.xCells; xBets = xs(betInds);
                ys = obj.yCells; yBets = ys(betInds);
                ns = obj.nCells; nBets = ns(betInds);
                ls = obj.lCells; lBets = ls(betInds);
                thets = obj.thetCells; thetBets = thets(betInds);
                
                %Get dynamics
                if length(obj.aCells) >= 1
                    xAlph = obj.xCells(i);
                    yAlph = obj.yCells(i);
                    lAlph = obj.lCells(i);
                    nAlph = obj.nCells(i);
                    thetAlph = obj.thetCells(i);
                    uAlph = obj.uCells(i,:)';
                end
                
                if compiled
                    [dUdx,dUdy,dUdthet] = mexCalcEnergyGradientsPeriodic(xBets,yBets,lBets,nBets,thetBets,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height);
                else
                    [dUdx,dUdy,dUdthet] = calcEnergyGradientsPeriodic(xBets,yBets,lBets,nBets,thetBets,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height);
                end
                
                gradXY = -[sum(dUdx)/2,sum(dUdy)/2];
                gradTheta = -sum(dUdthet)/2;
                
                v0 = obj.fCells(i)/(f0*fPar(i)); %The self-propulsion velocity of a non-interacting SPR
                
                drdt(i,:) = (v0*uAlph - (fT(:,:,i)\gradXY'))';
                dthetadt(i) = -gradTheta/fR(i);
            end
        end
        
        function obj = stepModel(obj,timeStep,f0,compiled)
            %Increases the time by one step, updating cell positions based on their current velocities.
            
            %Calculate distance matrix and threshold if needed (e.g. for first time point)
            if isempty(obj.distThresh)
                obj = obj.calcDistThresh();
            end
            if isempty(obj.cellDists)
                obj = obj.calcDistMat();
            end
                         
            %Apply midpoint method to simulate movement dynamics
            [drdtk1,dthetadtk1] = obj.calcVelocities(f0,compiled);
            k1 = obj.moveCells(drdtk1,dthetadtk1,timeStep/2);
            [drdt,dthetadt] = k1.calcVelocities(f0,compiled);

            %Update position and angle of cells based on dynamic parameters and timestep size
            obj = obj.moveCells(drdt,dthetadt,timeStep);
            
            %Calculate distance matrix. Allows elimination of small-intensity interactions at later time points.
            obj = obj.calcDistThresh();
            obj = obj.calcDistMat();
        end
        
        function obj = moveCells(obj,drdt,dthetadt,timeStep)
            %Updates the position of the cell based on current velocity
            obj.xCells = obj.xCells + timeStep*drdt(:,1);
            obj.yCells = obj.yCells + timeStep*drdt(:,2);
            
            %Apply periodic boundary conditions
            obj.xCells = mod(obj.xCells,obj.xWidth);
            obj.yCells = mod(obj.yCells,obj.yHeight);
                        
            %Updates the angles of the cell
            obj.thetCells = obj.thetCells + dthetadt*timeStep;
            obj.thetCells = mod(obj.thetCells + pi,2*pi) - pi;
            
            obj.uCells = [cos(obj.thetCells),sin(obj.thetCells)];
        end
        
        function outImg = drawField(obj)
            %Draws the current state of the model - location and angles of all rods in model. Coloured rods are motile cells.
            %Note that this is only a 2D projection for debugging. For proper imaging of the 3D system, use the paraview scripts.
            Upsample = 10; %Extent to which the 'design' image should be interpolated to create smoother graphics.
            Downsample = 2.5; %Extent to which the final image should be scaled down to save space.
            
            %Draw barrier (as an image) and break into separate rgb
            %channels
            background = imresize(ones(ceil(obj.yHeight),ceil(obj.xWidth)),Upsample); %We'll do the smoothing later
            imgr = imgaussfilt(double(background),Upsample*obj.lam);
            imgg = imgaussfilt(double(background),Upsample*obj.lam);
            imgb = imgaussfilt(double(background),Upsample*obj.lam);
            
            %Draw rods
            
            %Do separate paintjobs for each of the three colour
            %channels
            imgr = paintEllipse(imgr,obj.xCells,obj.yCells,obj.aCells/1.7,obj.lam*ones(size(obj.aCells))/1.1,-rad2deg(obj.thetCells),obj.cCells(:,1),1/Upsample);
            imgg = paintEllipse(imgg,obj.xCells,obj.yCells,obj.aCells/1.7,obj.lam*ones(size(obj.aCells))/1.1,-rad2deg(obj.thetCells),obj.cCells(:,2),1/Upsample);
            imgb = paintEllipse(imgb,obj.xCells,obj.yCells,obj.aCells/1.7,obj.lam*ones(size(obj.aCells))/1.1,-rad2deg(obj.thetCells),obj.cCells(:,3),1/Upsample);
            
            %Smooth to make nicer looking
            imgr = imgaussfilt(imgr,Upsample*obj.lam/4);
            imgg = imgaussfilt(imgg,Upsample*obj.lam/4);
            imgb = imgaussfilt(imgb,Upsample*obj.lam/4);
            
            outImg = cat(3,imgr,imgg,imgb);            
            outImg = imresize(outImg,1/Downsample);
        end
    end
end