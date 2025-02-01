classdef GIF_report_window < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        PleaseselectthemeshandLabel  matlab.ui.control.Label
        MeshStatisticsPanel          matlab.ui.container.Panel
        VerticesLabel                matlab.ui.control.Label
        TrianglesLabel               matlab.ui.control.Label
        BoundaryVerticesLabel        matlab.ui.control.Label
        MetaVerticesLabel            matlab.ui.control.Label
        verticesLabel                matlab.ui.control.Label
        trianglesLabel               matlab.ui.control.Label
        boundary_verticesLabel       matlab.ui.control.Label
        meta_verticesLabel           matlab.ui.control.Label
        ResultInformationPanel       matlab.ui.container.Panel
        E_SDAverageLabel             matlab.ui.control.Label
        E_SDAverageOver99oftheTrianglesWiththeLowestDistortionLabel  matlab.ui.control.Label
        kAverageLabel                matlab.ui.control.Label
        SegmentSizeLabel             matlab.ui.control.Label
        BoundaryMetaTrianglesLabel   matlab.ui.control.Label
        InternalMetaTrianglesLabel   matlab.ui.control.Label
        HastheFoldoversFixingMethodBeenAppliedLabel  matlab.ui.control.Label
        IstheFinalResultGloballyInjectiveLabel  matlab.ui.control.Label
        FoldoversinFinalResultLabel  matlab.ui.control.Label
        TotalTimesecondsLabel        matlab.ui.control.Label
        FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel  matlab.ui.control.Label
        totaltimeLabel               matlab.ui.control.Label
        segment_sizeLabel            matlab.ui.control.Label
        bmtLabel                     matlab.ui.control.Label
        imtLabel                     matlab.ui.control.Label
        eaLabel                      matlab.ui.control.Label
        ea99Label                    matlab.ui.control.Label
        fbLabel                      matlab.ui.control.Label
        keLabel                      matlab.ui.control.Label
        giLabel                      matlab.ui.control.Label
        ffLabel                      matlab.ui.control.Label
        hffmbaLabel                  matlab.ui.control.Label
        ScaffoldTriangulationsLabel  matlab.ui.control.Label
        InternalIterationsLabel      matlab.ui.control.Label
        sctLabel                     matlab.ui.control.Label
        iniLabel                     matlab.ui.control.Label
        VisualizeMeshinMatlabButton  matlab.ui.control.Button
        ExportResulttoobjFileButton  matlab.ui.control.Button
        VisualizeUVLayoutinMatlabButton  matlab.ui.control.Button
        ExportResulttomatFileButton  matlab.ui.control.Button
        SaveMeshStatisticsandResultInformationButton  matlab.ui.control.Button
    end

    
    properties (Access = public)
        mesh % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mesh)
            app.mesh = mesh;
            app.verticesLabel.Text = num2str(mesh.parametersMatrix(1,1));
            app.trianglesLabel.Text = num2str(mesh.parametersMatrix(2,1));
            app.boundary_verticesLabel.Text = num2str(mesh.parametersMatrix(3,1));
            app.meta_verticesLabel.Text = num2str(mesh.parametersMatrix(4,1));
            app.totaltimeLabel.Text = num2str(mesh.parametersMatrix(5,1));
            app.segment_sizeLabel.Text = num2str(mesh.parametersMatrix(6,1));
            app.bmtLabel.Text = num2str(mesh.parametersMatrix(7,1));
            app.imtLabel.Text = num2str(mesh.parametersMatrix(8,1));
            app.eaLabel.Text = num2str(mesh.parametersMatrix(9,1),'%.12f');
            app.ea99Label.Text = num2str(mesh.parametersMatrix(10,1),'%.12f');
            app.keLabel.Text = num2str(mesh.parametersMatrix(11,1),'%.7f');
            if mesh.parametersMatrix(12,1) == 1
                app.giLabel.Text = 'Yes';
            else 
                app.giLabel.Text = 'No';
            end
            app.ffLabel.Text = num2str(mesh.parametersMatrix(13,1));
            if mesh.parametersMatrix(14,1) == 1
                app.hffmbaLabel.Text = 'Yes';
                app.FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel.Visible = 'on';
                app.fbLabel.Visible = 'on';
            else 
                app.hffmbaLabel.Text = 'No';
                app.FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel.Visible = 'off';
                app.fbLabel.Visible = 'off';
            end
            app.fbLabel.Text = num2str(mesh.parametersMatrix(15,1));
            app.sctLabel.Text = num2str(mesh.parametersMatrix(16,1));
            app.iniLabel.Text = num2str(mesh.parametersMatrix(17,1));
        end

        % Button pushed function: VisualizeMeshinMatlabButton
        function VisualizeMeshinMatlabButtonPushed(app, event)
            visualizeMesh(app.mesh);
        end

        % Button pushed function: VisualizeUVLayoutinMatlabButton
        function VisualizeUVLayoutinMatlabButtonPushed(app, event)
            visualizeUVLayout(app.mesh)
        end

        % Button pushed function: ExportResulttoobjFileButton
        function ExportResulttoobjFileButtonPushed(app, event)
            choice = questdlg(['Saving the result as an .obj file might lead to errors in extreme cases, due to the conversion of double to ASCII.' ...
                'Do you want to continue?'], 'Warning','Yes','No','Yes');
            if strcmp(choice,'No') == 1
                return;
            end
            [name,path] = uiputfile(['.obj']);
            if name == 0
                return
            end
            fileID = fopen([path name],'w+t');
            saveObj( path, name, app.mesh );
        end

        % Button pushed function: ExportResulttomatFileButton
        function ExportResulttomatFileButtonPushed(app, event)
            [name,path] = uiputfile(['.mat']);
            if name == 0
                return
            end
            V = app.mesh.V;
            F = app.mesh.F1;
            UVs = app.mesh.uvs;
            save([path name],"V","F", "UVs");
        end

        % Button pushed function: 
        % SaveMeshStatisticsandResultInformationButton
        function SaveMeshStatisticsandResultInformationButtonPushed(app, event)
            [name,path] = uiputfile(['.txt']);
            if name == 0
                return
            end
            fileID = fopen([path name],'w+t');
            fprintf(fileID,'*******Globally Injective Flattening via a Reduced Harmonic Subspace - log******\n');
            fprintf(fileID,'Model name: %s\n', app.mesh.meshName);
            fprintf(fileID,'#Vertices: %i\n',app.mesh.parametersMatrix(1,1));
            fprintf(fileID,'#Faces: %i\n',app.mesh.parametersMatrix(2,1));
            fprintf(fileID,'#Boundary vertices: %i\n',app.mesh.parametersMatrix(3,1));
            fprintf(fileID,'#Meta vertices: %i\n',app.mesh.parametersMatrix(4,1));
            fprintf(fileID,'Total time %i\n',app.mesh.parametersMatrix(5,1));
            fprintf(fileID,'Segment size %i\n',app.mesh.parametersMatrix(6,1));
            fprintf(fileID,'#Boundary meta triangles: %i\n',app.mesh.parametersMatrix(7,1));
            fprintf(fileID,'#Internal meta triangles: %i\n',app.mesh.parametersMatrix(8,1));
            fprintf(fileID,'E_SD average: %f\n',app.mesh.parametersMatrix(9,1));
            fprintf(fileID,'E_SD average over 99 percents of the triangles with the lowest distortion: %f\n',app.mesh.parametersMatrix(10,1));
            fprintf(fileID,'k Average: %i\n',app.mesh.parametersMatrix(11,1));
            fprintf(fileID,'#Scaffold triangulations: %i\n',app.mesh.parametersMatrix(16,1));
            fprintf(fileID,'#Internal iterations: %i\n',app.mesh.parametersMatrix(17,1));
            if app.mesh.parametersMatrix(12,1) == 1
                fprintf(fileID,'The final result is globally injective\n');
            else 
                fprintf(fileID,'The final result is not globally injective\n');
            end
            fprintf(fileID,'#Foldovers in final result: %i\n',app.mesh.parametersMatrix(13,1));
            if app.mesh.parametersMatrix(14,1) == 1
                fprintf(fileID,'The foldovers fixing method has been applied\n');
                fprintf(fileID,'#Foldovers before the foldovers fixing method was applied: %i\n',app.mesh.parametersMatrix(15,1));
            else 
                fprintf(fileID,'The foldovers fixing method has not been applied\n');
            end
            fclose('all');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1238 759];
            app.UIFigure.Name = 'UI Figure';

            % Create PleaseselectthemeshandLabel
            app.PleaseselectthemeshandLabel = uilabel(app.UIFigure);
            app.PleaseselectthemeshandLabel.HorizontalAlignment = 'center';
            app.PleaseselectthemeshandLabel.FontSize = 24;
            app.PleaseselectthemeshandLabel.Position = [200 656 793 30];
            app.PleaseselectthemeshandLabel.Text = 'Globally Injective Flattening via a Reduced Harmonic Subspace - Report';

            % Create MeshStatisticsPanel
            app.MeshStatisticsPanel = uipanel(app.UIFigure);
            app.MeshStatisticsPanel.Title = 'Mesh Statistics:';
            app.MeshStatisticsPanel.FontSize = 19;
            app.MeshStatisticsPanel.Position = [57 401 468 233];

            % Create VerticesLabel
            app.VerticesLabel = uilabel(app.MeshStatisticsPanel);
            app.VerticesLabel.FontSize = 19;
            app.VerticesLabel.Position = [21 158 89 23];
            app.VerticesLabel.Text = '#Vertices:';

            % Create TrianglesLabel
            app.TrianglesLabel = uilabel(app.MeshStatisticsPanel);
            app.TrianglesLabel.FontSize = 19;
            app.TrianglesLabel.Position = [21 119 99 23];
            app.TrianglesLabel.Text = '#Triangles:';

            % Create BoundaryVerticesLabel
            app.BoundaryVerticesLabel = uilabel(app.MeshStatisticsPanel);
            app.BoundaryVerticesLabel.FontSize = 19;
            app.BoundaryVerticesLabel.Position = [22 76 176 23];
            app.BoundaryVerticesLabel.Text = '#Boundary Vertices:';

            % Create MetaVerticesLabel
            app.MetaVerticesLabel = uilabel(app.MeshStatisticsPanel);
            app.MetaVerticesLabel.FontSize = 19;
            app.MetaVerticesLabel.Position = [21 32 136 23];
            app.MetaVerticesLabel.Text = '#Meta Vertices:';

            % Create verticesLabel
            app.verticesLabel = uilabel(app.MeshStatisticsPanel);
            app.verticesLabel.FontSize = 19;
            app.verticesLabel.Position = [293 158 161 23];
            app.verticesLabel.Text = 'vertices';

            % Create trianglesLabel
            app.trianglesLabel = uilabel(app.MeshStatisticsPanel);
            app.trianglesLabel.FontSize = 19;
            app.trianglesLabel.Position = [293 119 161 23];
            app.trianglesLabel.Text = 'triangles';

            % Create boundary_verticesLabel
            app.boundary_verticesLabel = uilabel(app.MeshStatisticsPanel);
            app.boundary_verticesLabel.FontSize = 19;
            app.boundary_verticesLabel.Position = [293 76 161 23];
            app.boundary_verticesLabel.Text = 'boundary_vertices';

            % Create meta_verticesLabel
            app.meta_verticesLabel = uilabel(app.MeshStatisticsPanel);
            app.meta_verticesLabel.FontSize = 19;
            app.meta_verticesLabel.Position = [293 32 161 23];
            app.meta_verticesLabel.Text = 'meta_vertices';

            % Create ResultInformationPanel
            app.ResultInformationPanel = uipanel(app.UIFigure);
            app.ResultInformationPanel.Title = 'Result Information';
            app.ResultInformationPanel.FontSize = 19;
            app.ResultInformationPanel.Position = [600 22 610 612];

            % Create E_SDAverageLabel
            app.E_SDAverageLabel = uilabel(app.ResultInformationPanel);
            app.E_SDAverageLabel.FontSize = 19;
            app.E_SDAverageLabel.Position = [19 379 141 23];
            app.E_SDAverageLabel.Text = 'E_SD Average :';

            % Create E_SDAverageOver99oftheTrianglesWiththeLowestDistortionLabel
            app.E_SDAverageOver99oftheTrianglesWiththeLowestDistortionLabel = uilabel(app.ResultInformationPanel);
            app.E_SDAverageOver99oftheTrianglesWiththeLowestDistortionLabel.FontSize = 19;
            app.E_SDAverageOver99oftheTrianglesWiththeLowestDistortionLabel.Position = [19 321 318 44];
            app.E_SDAverageOver99oftheTrianglesWiththeLowestDistortionLabel.Text = {'E_SD Average Over 99% of the'; 'Triangles With the Lowest Distortion:'};

            % Create kAverageLabel
            app.kAverageLabel = uilabel(app.ResultInformationPanel);
            app.kAverageLabel.FontSize = 19;
            app.kAverageLabel.Position = [19 286 96 23];
            app.kAverageLabel.Text = 'k Average:';

            % Create SegmentSizeLabel
            app.SegmentSizeLabel = uilabel(app.ResultInformationPanel);
            app.SegmentSizeLabel.FontSize = 19;
            app.SegmentSizeLabel.Position = [19 498 129 23];
            app.SegmentSizeLabel.Text = 'Segment Size:';

            % Create BoundaryMetaTrianglesLabel
            app.BoundaryMetaTrianglesLabel = uilabel(app.ResultInformationPanel);
            app.BoundaryMetaTrianglesLabel.FontSize = 19;
            app.BoundaryMetaTrianglesLabel.Position = [19 459 235 23];
            app.BoundaryMetaTrianglesLabel.Text = '#Boundary Meta Triangles:';

            % Create InternalMetaTrianglesLabel
            app.InternalMetaTrianglesLabel = uilabel(app.ResultInformationPanel);
            app.InternalMetaTrianglesLabel.FontSize = 19;
            app.InternalMetaTrianglesLabel.Position = [19 416 216 23];
            app.InternalMetaTrianglesLabel.Text = '#Internal Meta Triangles:';

            % Create HastheFoldoversFixingMethodBeenAppliedLabel
            app.HastheFoldoversFixingMethodBeenAppliedLabel = uilabel(app.ResultInformationPanel);
            app.HastheFoldoversFixingMethodBeenAppliedLabel.FontSize = 19;
            app.HastheFoldoversFixingMethodBeenAppliedLabel.Position = [19 76 217 44];
            app.HastheFoldoversFixingMethodBeenAppliedLabel.Text = {'Has the Foldovers Fixing'; 'Method Been Applied?'};

            % Create IstheFinalResultGloballyInjectiveLabel
            app.IstheFinalResultGloballyInjectiveLabel = uilabel(app.ResultInformationPanel);
            app.IstheFinalResultGloballyInjectiveLabel.FontSize = 19;
            app.IstheFinalResultGloballyInjectiveLabel.Position = [19 179 321 23];
            app.IstheFinalResultGloballyInjectiveLabel.Text = 'Is the Final Result Globally Injective?';

            % Create FoldoversinFinalResultLabel
            app.FoldoversinFinalResultLabel = uilabel(app.ResultInformationPanel);
            app.FoldoversinFinalResultLabel.FontSize = 19;
            app.FoldoversinFinalResultLabel.Position = [19 140 232 23];
            app.FoldoversinFinalResultLabel.Text = '#Foldovers in Final Result:';

            % Create TotalTimesecondsLabel
            app.TotalTimesecondsLabel = uilabel(app.ResultInformationPanel);
            app.TotalTimesecondsLabel.FontSize = 19;
            app.TotalTimesecondsLabel.Position = [19 537 186 23];
            app.TotalTimesecondsLabel.Text = 'Total Time (seconds):';

            % Create FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel
            app.FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel = uilabel(app.ResultInformationPanel);
            app.FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel.FontSize = 19;
            app.FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel.Position = [19 15 284 44];
            app.FoldoversBeforetheFoldoversFixingMethodWasAppliedLabel.Text = {'#Foldovers Before the Foldovers'; ' Fixing Method Was Applied:'};

            % Create totaltimeLabel
            app.totaltimeLabel = uilabel(app.ResultInformationPanel);
            app.totaltimeLabel.FontSize = 19;
            app.totaltimeLabel.Position = [375 537 153 23];
            app.totaltimeLabel.Text = 'totaltime';

            % Create segment_sizeLabel
            app.segment_sizeLabel = uilabel(app.ResultInformationPanel);
            app.segment_sizeLabel.FontSize = 19;
            app.segment_sizeLabel.Position = [376 498 202 23];
            app.segment_sizeLabel.Text = 'segment_size';

            % Create bmtLabel
            app.bmtLabel = uilabel(app.ResultInformationPanel);
            app.bmtLabel.FontSize = 19;
            app.bmtLabel.Position = [376 459 131 23];
            app.bmtLabel.Text = 'bmt';

            % Create imtLabel
            app.imtLabel = uilabel(app.ResultInformationPanel);
            app.imtLabel.FontSize = 19;
            app.imtLabel.Position = [375 420 132 23];
            app.imtLabel.Text = 'imt';

            % Create eaLabel
            app.eaLabel = uilabel(app.ResultInformationPanel);
            app.eaLabel.FontSize = 19;
            app.eaLabel.Position = [376 381 214 23];
            app.eaLabel.Text = 'ea';

            % Create ea99Label
            app.ea99Label = uilabel(app.ResultInformationPanel);
            app.ea99Label.FontSize = 19;
            app.ea99Label.Position = [376 332 214 23];
            app.ea99Label.Text = 'ea99';

            % Create fbLabel
            app.fbLabel = uilabel(app.ResultInformationPanel);
            app.fbLabel.FontSize = 19;
            app.fbLabel.Position = [375 30 77 23];
            app.fbLabel.Text = 'fb';

            % Create keLabel
            app.keLabel = uilabel(app.ResultInformationPanel);
            app.keLabel.FontSize = 19;
            app.keLabel.Position = [376 284 152 23];
            app.keLabel.Text = 'ke';

            % Create giLabel
            app.giLabel = uilabel(app.ResultInformationPanel);
            app.giLabel.FontSize = 19;
            app.giLabel.Position = [375 179 88 23];
            app.giLabel.Text = 'gi';

            % Create ffLabel
            app.ffLabel = uilabel(app.ResultInformationPanel);
            app.ffLabel.FontSize = 19;
            app.ffLabel.Position = [376 140 63 23];
            app.ffLabel.Text = 'ff';

            % Create hffmbaLabel
            app.hffmbaLabel = uilabel(app.ResultInformationPanel);
            app.hffmbaLabel.FontSize = 19;
            app.hffmbaLabel.Position = [375 87 77 23];
            app.hffmbaLabel.Text = 'hffmba';

            % Create ScaffoldTriangulationsLabel
            app.ScaffoldTriangulationsLabel = uilabel(app.ResultInformationPanel);
            app.ScaffoldTriangulationsLabel.FontSize = 19;
            app.ScaffoldTriangulationsLabel.Position = [19 251 215 23];
            app.ScaffoldTriangulationsLabel.Text = '#Scaffold Triangulations:';

            % Create InternalIterationsLabel
            app.InternalIterationsLabel = uilabel(app.ResultInformationPanel);
            app.InternalIterationsLabel.FontSize = 19;
            app.InternalIterationsLabel.Position = [19 215 168 23];
            app.InternalIterationsLabel.Text = '#Internal Iterations:';

            % Create sctLabel
            app.sctLabel = uilabel(app.ResultInformationPanel);
            app.sctLabel.FontSize = 19;
            app.sctLabel.Position = [376 251 131 23];
            app.sctLabel.Text = 'sct';

            % Create iniLabel
            app.iniLabel = uilabel(app.ResultInformationPanel);
            app.iniLabel.FontSize = 19;
            app.iniLabel.Position = [375 215 88 23];
            app.iniLabel.Text = 'ini';

            % Create VisualizeMeshinMatlabButton
            app.VisualizeMeshinMatlabButton = uibutton(app.UIFigure, 'push');
            app.VisualizeMeshinMatlabButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizeMeshinMatlabButtonPushed, true);
            app.VisualizeMeshinMatlabButton.FontSize = 19;
            app.VisualizeMeshinMatlabButton.Position = [57 323 264 53];
            app.VisualizeMeshinMatlabButton.Text = 'Visualize Mesh in Matlab';

            % Create ExportResulttoobjFileButton
            app.ExportResulttoobjFileButton = uibutton(app.UIFigure, 'push');
            app.ExportResulttoobjFileButton.ButtonPushedFcn = createCallbackFcn(app, @ExportResulttoobjFileButtonPushed, true);
            app.ExportResulttoobjFileButton.FontSize = 19;
            app.ExportResulttoobjFileButton.Position = [57 171 264 53];
            app.ExportResulttoobjFileButton.Text = 'Export Result to .obj File';

            % Create VisualizeUVLayoutinMatlabButton
            app.VisualizeUVLayoutinMatlabButton = uibutton(app.UIFigure, 'push');
            app.VisualizeUVLayoutinMatlabButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizeUVLayoutinMatlabButtonPushed, true);
            app.VisualizeUVLayoutinMatlabButton.FontSize = 19;
            app.VisualizeUVLayoutinMatlabButton.Position = [57 245 264 53];
            app.VisualizeUVLayoutinMatlabButton.Text = 'Visualize UV Layout in Matlab';

            % Create ExportResulttomatFileButton
            app.ExportResulttomatFileButton = uibutton(app.UIFigure, 'push');
            app.ExportResulttomatFileButton.ButtonPushedFcn = createCallbackFcn(app, @ExportResulttomatFileButtonPushed, true);
            app.ExportResulttomatFileButton.FontSize = 19;
            app.ExportResulttomatFileButton.Position = [57 98 264 53];
            app.ExportResulttomatFileButton.Text = 'Export Result to .mat File';

            % Create SaveMeshStatisticsandResultInformationButton
            app.SaveMeshStatisticsandResultInformationButton = uibutton(app.UIFigure, 'push');
            app.SaveMeshStatisticsandResultInformationButton.ButtonPushedFcn = createCallbackFcn(app, @SaveMeshStatisticsandResultInformationButtonPushed, true);
            app.SaveMeshStatisticsandResultInformationButton.FontSize = 19;
            app.SaveMeshStatisticsandResultInformationButton.Position = [57 22 264 53];
            app.SaveMeshStatisticsandResultInformationButton.Text = {'Save Mesh Statistics'; 'and Result Information'};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GIF_report_window(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end