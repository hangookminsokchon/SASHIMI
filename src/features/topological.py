"""
Description: Compute a suite of topological features 

Dependencies: Refer ReadMe.

List of topological features
- Cubical Complex
    - Tumor-Immune, Tumor-Stromal, Immune-Stroaml pairs

- Witness Complex
    - Tumor-Immune, Tumor-Stromal, Immune-Stroaml pairs
"""
#-------------------------------------------------------------------------------
import pandas as pd                             # Dataframe manipulation
import numpy as np                              # Matrix operations
import gudhi                                    # Topological feature computation ()
import matplotlib.pyplot as plt                 # Data visualization
from matplotlib.colors import ListedColormap    # =
from scipy.stats import gaussian_kde            # Gaussian KDE for type assign
from scipy import ndimage                       # 
from scipy.ndimage import distance_transform_bf # 
#-------------------------------------------------------------------------------

class TopoFeature:
    # Constructor
    def __init__(self):
        # Cell types
        self.ordered_celltype = ['immune', 'stromal', 'tumor']
        # Colors by cell types
        self.ordered_colortype = ['black', 'red', 'green']
        # Cell type pairs
        self.cell_pairs = [
            ('tumor', 'immune'),
            ('tumor', 'stromal'),
            ('immune', 'stromal')
        ]
        
    # Read in raw data, assigns the cell types
    def read_img(self, img_dir):
        """
        Read the raw image and convert them into labeled point pattern data
        
        Parameters:
        -----------
        img_dir : str
            Path of input image
            
        Returns:
        -----------
        points : np.array
            Point pattern data of input image
        labels : np.array
            Point pattern data's corresponding cell type labels
        """
        df = pd.read_csv(img_dir)
        df.columns = ['x', 'y', 'type']
    
        points = df.iloc[:, 0:2].values.astype(np.float64)
        labels = df['type'].astype(str)
    
        return [points, labels]       
    
    
    # Gaussian KDE type classifier for cubical complex based on SEDT3
    def compute_kde(self, points, labels):
        """
        Assign cell types to empty grid based on densities computed by Gaussian KDE
        
        Parameters:
        -----------
        points : np.array
            Point pattern data of input image
        labels : np.array
            Point pattern data's corresponding cell type labels
            
        Returns:
        --------
        class_map : np.array
            Grid with assinged cell types
        """
        label_to_index = {lbl: i for i, lbl in enumerate(self.ordered_celltype)}
        n_labels = len(self.ordered_celltype)
        
        # Evaluate gaussian kde on 100 x 100 grid
        
        # Generate 100 x 100 empty grid
        grid_res = 100
        x, y = np.linspace(0, 1, grid_res), np.linspace(0, 1, grid_res)
        xx, yy = np.meshgrid(x, y)
        grid_coords = np.vstack([xx.ravel(), yy.ravel()])
        
        grid_coords = np.asarray(grid_coords, dtype = np.float64)
        
        # KDE
        densities = np.zeros((n_labels, grid_coords.shape[1]))
        
        for i, label in enumerate(self.ordered_celltype):
            pts = points[labels == label]
            if len(pts) > 0:  # Check if there are points for this label
                kde = gaussian_kde(pts.T, bw_method = 0.5)
                densities[i] = kde(grid_coords)
            
        # Assigning class by density in each grid
        class_map = np.argmax(densities, axis = 0).reshape((grid_res, grid_res))
        
        return class_map
    
    
    def compute_cubical_complex_pair(self, class_map, cell_type1, cell_type2):
        """
        Compute cubical complex for a specific cell type pair
        
        Parameters:
        -----------
        class_map : np.array
            Class map from KDE
        cell_type1, cell_type2 : str
            Cell types to compare
            
         Returns:
        --------
        cubical_complex : np.array
            
        diag : list
            Persistence diagram from GUDHI
        """
        
        # Map cell types to indices
        type_to_idx = {'immune': 0, 'stromal': 1, 'tumor': 2}
        idx1 = type_to_idx[cell_type1]
        idx2 = type_to_idx[cell_type2]
        
        # SEDT3 computation
        class_grid = pd.DataFrame(class_map)
        
        # Create binary grids for each cell type
        grid_type1 = (class_grid == idx1).astype(np.int32)
        grid_type2 = (class_grid == idx2).astype(np.int32)
        
        # Compute distance transforms
        dist_type1 = distance_transform_bf(grid_type1.values, metric='euclidean').astype(np.float64)
        dist_type2 = distance_transform_bf(1 - grid_type2.values, metric='euclidean').astype(np.float64)
        
        # SEDT3
        sedt3 = dist_type2 - dist_type1
        
        # Set regions of other cell types to infinity
        other_idx = [i for i in range(3) if i not in [idx1, idx2]]
        for idx in other_idx:
            sedt3[class_grid == idx] = np.inf
        
        dims = list(sedt3.shape)
        flat_vals = sedt3.flatten()
        
        cubical_complex = gudhi.CubicalComplex(dimensions = dims, top_dimensional_cells = flat_vals)
        diag = cubical_complex.persistence()
        
        return [cubical_complex, diag]
    
    
    def compute_witness_complex_pair(self, points, labels, cell_type1, cell_type2):
        """
        Compute witness complex for a specific cell type pair
        
        Parameters:
        -----------
        points : np.array
            Point coordinates
        labels : np.array
            Cell type labels
        cell_type1, cell_type2 : str
            Cell types to compare (landmarks and witnesses)

        Returns:
        --------
         diag : list
            Persistence diagram from GUDHI
        
        """
        
        points_type1 = points[labels == cell_type1]
        points_type2 = points[labels == cell_type2]
        
        # Check if both cell types have enough points
        if len(points_type1) < 3 or len(points_type2) < 3:
            # Return empty diagram if not enough points
            return [None, None, []]
        
        # Determine number of points to sample
        if max(len(points_type1), len(points_type2)) <= 2000:
            nb_points = int(max(len(points_type1), len(points_type2)) / 10)
        else:
            nb_points = 200
        
        # Ensure we don't sample more points than available
        nb_points = min(nb_points, len(points_type1), len(points_type2))
        nb_points = max(nb_points, 3)  # At least 3 points
        
        max_alpha_square = 0.1
        limit_dimension = 2
        
        # Sample landmarks and witnesses
        landmarks = gudhi.subsampling.pick_n_random_points(points=points_type1, nb_points=nb_points)
        witnesses = gudhi.subsampling.pick_n_random_points(points=points_type2, nb_points=nb_points)
        
        witness_complex = gudhi.EuclideanWitnessComplex(witnesses=witnesses, landmarks=landmarks)
        simplex_tree = witness_complex.create_simplex_tree(max_alpha_square=max_alpha_square, 
                                                           limit_dimension=limit_dimension)
        diag = simplex_tree.persistence()
        
        return [witness_complex, simplex_tree, diag]
    
    
    def extract_persistence_stats(self, diag, prefix='', exclude_infinite=True, max_finite_value=None):
        """
        Extract summary statistics from persistence diagram
        
        Parameters:
        -----------
        diag : list
            Persistence diagram from GUDHI
        prefix : str
            Prefix for column names
        exclude_infinite : bool
            Whether to exclude points with infinite death time
        max_finite_value : float
            Replace infinite values with this value (if exclude_infinite=False)
            
        Returns:
        --------
        stats: dict
            Dictionary of statistics
        """
        
        # Handle empty diagram
        if not diag or len(diag) == 0:
            return {
                f'{prefix}_birth_min': np.nan,
                f'{prefix}_birth_max': np.nan,
                f'{prefix}_birth_mean': np.nan,
                f'{prefix}_birth_std': np.nan,
                f'{prefix}_death_min': np.nan,
                f'{prefix}_death_max': np.nan,
                f'{prefix}_death_mean': np.nan,
                f'{prefix}_death_std': np.nan,
                f'{prefix}_lifetime_min': np.nan,
                f'{prefix}_lifetime_max': np.nan,
                f'{prefix}_lifetime_mean': np.nan,
                f'{prefix}_lifetime_std': np.nan,
                f'{prefix}_n_features': 0
            }
        
        # Extract birth and death times
        births = []
        deaths = []
        
        for interval in diag:
            dim, (birth, death) = interval
            
            # Handle infinite values
            if np.isinf(death):
                if exclude_infinite:
                    continue
                elif max_finite_value is not None:
                    death = max_finite_value
            
            births.append(birth)
            deaths.append(death)
        
        # Convert to numpy arrays
        births = np.array(births) if births else np.array([])
        deaths = np.array(deaths) if deaths else np.array([])
        
        # Calculate statistics
        if len(births) > 0:
            lifetimes = deaths - births
            stats = {
                f'{prefix}_birth_min': np.min(births),
                f'{prefix}_birth_max': np.max(births),
                f'{prefix}_birth_mean': np.mean(births),
                f'{prefix}_birth_std': np.std(births),
                f'{prefix}_death_min': np.min(deaths),
                f'{prefix}_death_max': np.max(deaths),
                f'{prefix}_death_mean': np.mean(deaths),
                f'{prefix}_death_std': np.std(deaths),
                f'{prefix}_lifetime_min': np.min(lifetimes),
                f'{prefix}_lifetime_max': np.max(lifetimes),
                f'{prefix}_lifetime_mean': np.mean(lifetimes),
                f'{prefix}_lifetime_std': np.std(lifetimes),
                f'{prefix}_n_features': len(births)
            }
        else:
            stats = {
                f'{prefix}_birth_min': np.nan,
                f'{prefix}_birth_max': np.nan,
                f'{prefix}_birth_mean': np.nan,
                f'{prefix}_birth_std': np.nan,
                f'{prefix}_death_min': np.nan,
                f'{prefix}_death_max': np.nan,
                f'{prefix}_death_mean': np.nan,
                f'{prefix}_death_std': np.nan,
                f'{prefix}_lifetime_min': np.nan,
                f'{prefix}_lifetime_max': np.nan,
                f'{prefix}_lifetime_mean': np.nan,
                f'{prefix}_lifetime_std': np.nan,
                f'{prefix}_n_features': 0
            }
        
        return stats
    
    
    def compute_all_features(self, img_dir, exclude_infinite=True):
        """
        Compute all topological features for all cell type pairs and return as 1-row DataFrame
        
        Parameters:
        -----------
        img_dir : str
            Path to image CSV file
        exclude_infinite : bool
            Whether to exclude infinite persistence values
            
        Returns:
        --------
        df_features : pd.DataFrame
            1-row DataFrame with all features
        """
        
        # Read image
        points, labels = self.read_img(img_dir)
        
        # Compute KDE once for all cubical complexes
        class_map = self.compute_kde(points, labels)
        
        # Initialize stats dictionary
        all_stats = {}
        
        # Compute features for each cell type pair
        for cell_type1, cell_type2 in self.cell_pairs:
            pair_name = f'{cell_type1}_{cell_type2}'
            
            # Cubical complex features
            try:
                _, cubical_diag = self.compute_cubical_complex_pair(class_map, cell_type1, cell_type2)
                cubical_stats = self.extract_persistence_stats(
                    cubical_diag, 
                    prefix=f'cubical_{pair_name}',
                    exclude_infinite=exclude_infinite
                )
                all_stats.update(cubical_stats)
            except Exception as e:
                print(f"Error computing cubical complex for {pair_name}: {e}")
                # Add NaN values for this pair
                cubical_stats = self.extract_persistence_stats(
                    [], 
                    prefix=f'cubical_{pair_name}',
                    exclude_infinite=exclude_infinite
                )
                all_stats.update(cubical_stats)
            
            # Witness complex features
            try:
                _, _, witness_diag = self.compute_witness_complex_pair(points, labels, cell_type1, cell_type2)
                witness_stats = self.extract_persistence_stats(
                    witness_diag, 
                    prefix=f'witness_{pair_name}',
                    exclude_infinite=exclude_infinite
                )
                all_stats.update(witness_stats)
            except Exception as e:
                print(f"Error computing witness complex for {pair_name}: {e}")
                # Add NaN values for this pair
                witness_stats = self.extract_persistence_stats(
                    [], 
                    prefix=f'witness_{pair_name}',
                    exclude_infinite=exclude_infinite
                )
                all_stats.update(witness_stats)
        
        # Convert to 1-row DataFrame
        df_features = pd.DataFrame([all_stats])
        
        return df_features
    
    
    def compute_features_by_dimension(self, img_dir, dimensions=[0, 1], exclude_infinite=True):
        """
        Compute features separated by homology dimension for all cell type pairs
        
        Parameters:
        -----------
        img_dir : str
            Path to image CSV file
        dimensions : list
            List of dimensions to include (e.g., [0, 1] for H0 and H1)
        exclude_infinite : bool
            Whether to exclude infinite persistence values
            
        Returns:
        --------
        df_features : pd.DataFrame
            1-row DataFrame with features separated by dimension
        """
        
        # Read image
        points, labels = self.read_img(img_dir)
        
        # Compute KDE once for all cubical complexes
        class_map = self.compute_kde(points, labels)
        
        # Initialize stats dictionary
        all_stats = {}
        
        # Compute features for each cell type pair
        for cell_type1, cell_type2 in self.cell_pairs:
            pair_name = f'{cell_type1}_{cell_type2}'
            
            # Cubical complex features by dimension
            try:
                _, cubical_diag = self.compute_cubical_complex_pair(class_map, cell_type1, cell_type2)
                
                for dim in dimensions:
                    # Filter diagram by dimension
                    dim_diag = [(d, interval) for d, interval in cubical_diag if d == dim]
                    dim_stats = self.extract_persistence_stats(
                        dim_diag,
                        prefix=f'cubical_{pair_name}_h{dim}',
                        exclude_infinite=exclude_infinite
                    )
                    all_stats.update(dim_stats)
            except Exception as e:
                print(f"Error computing cubical complex for {pair_name}: {e}")
                for dim in dimensions:
                    dim_stats = self.extract_persistence_stats(
                        [],
                        prefix=f'cubical_{pair_name}_h{dim}',
                        exclude_infinite=exclude_infinite
                    )
                    all_stats.update(dim_stats)
            
            # Witness complex features by dimension
            try:
                _, _, witness_diag = self.compute_witness_complex_pair(points, labels, cell_type1, cell_type2)
                
                for dim in dimensions:
                    # Filter diagram by dimension
                    dim_diag = [(d, interval) for d, interval in witness_diag if d == dim]
                    dim_stats = self.extract_persistence_stats(
                        dim_diag,
                        prefix=f'witness_{pair_name}_h{dim}',
                        exclude_infinite=exclude_infinite
                    )
                    all_stats.update(dim_stats)
            except Exception as e:
                print(f"Error computing witness complex for {pair_name}: {e}")
                for dim in dimensions:
                    dim_stats = self.extract_persistence_stats(
                        [],
                        prefix=f'witness_{pair_name}_h{dim}',
                        exclude_infinite=exclude_infinite
                    )
                    all_stats.update(dim_stats)
        
        # Convert to 1-row DataFrame
        df_features = pd.DataFrame([all_stats])
        
        return df_features


