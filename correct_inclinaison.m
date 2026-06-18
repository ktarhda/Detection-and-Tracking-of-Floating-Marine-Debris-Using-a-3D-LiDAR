function ptCloud_corrige = correct_inclinaison(ptCloud_in, imu)

    % 1. Extraction et moyenne de l'accéléromètre
    if isstruct(imu) || istable(imu)
        accel_data = imu.AccelerometerReadings{:, :};
        ax = mean(accel_data(:,1));
        ay = mean(accel_data(:,2));
        az = mean(accel_data(:,3));
    else
        % Cas vecteur [ax, ay, az] direct (temps réel)
        ax = imu(1);
        ay = imu(2);
        az = imu(3);
    end

    % 2. Calcul Roll et Pitch (Yaw ignoré)
    roll  = atan2(ay, az);
    pitch = atan2(-ax, sqrt(ay^2 + az^2));

    Rx = [1,         0,          0;
          0,  cos(roll), -sin(roll);
          0,  sin(roll),  cos(roll)];

    Ry = [cos(pitch), 0, sin(pitch);
          0,          1,          0;
         -sin(pitch), 0, cos(pitch)];

    % Rotation inverse = correction de l'inclinaison
    R_correction = (Ry * Rx)';

    % 3. Translations JSON 
    % ousterFileReader applique déjà lidar_to_sensor en interne
    % → on ne fait QUE la correction IMU to  Sensor
    t_imu_sensor   = [0.006253, -0.011775, 0.007645] / 1000;
    t_lidar_sensor = [0, 0, 36.18] / 1000;

    % Décalage entre le centre LiDAR et l'IMU
    decalage_lidar_imu = t_imu_sensor - t_lidar_sensor;

    % 4. Trois transformations rigides
    % T1 : déplacer le nuage au centre de l'IMU
    tform_vers_imu = rigidtform3d(eye(3), -decalage_lidar_imu);

    % T2 : appliquer la rotation de correction
    tform_rotation = rigidtform3d(R_correction, [0, 0, 0]);

    % T3 : ramener au centre du LiDAR
    tform_retour = rigidtform3d(eye(3), decalage_lidar_imu);

    % 5. Application
    ptCloud_decale  = pctransform(ptCloud_in,    tform_vers_imu);
    ptCloud_tourne  = pctransform(ptCloud_decale, tform_rotation);
    ptCloud_corrige = pctransform(ptCloud_tourne, tform_retour);

end