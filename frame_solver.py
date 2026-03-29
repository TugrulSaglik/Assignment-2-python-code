"""
Module: frame_solver
Purpose: Handles the structural analysis pipeline utilizing the custom matrix library.
"""

from matrix_library import Matrix, SparseMatrix

class FrameModel:
    def __init__(self):
        self.materials = {}  
        self.nodes = {}      
        self.members = {}    
        self.supports = {}   
        self.loads = {}      
        
        self.equation_map = {}          
        self.num_equations = 0       
        
        self.k_local = {}   
        self.rot_matrices = {}         
        self.k_global_element = {}  
        
        self.K_global_struct = None         
        self.load_vector = []         
        self.displacements = []         
        self.member_forces = {}      

    def process_equations(self):
        """Assigns active degrees of freedom to equations."""
        self.equation_map = {n_id: [0, 0, 0] for n_id in self.nodes.keys()}
        
        for sup in self.supports.values():
            n_id = sup['node_id']
            self.equation_map[n_id] = [sup['RX'], sup['RY'], sup['RZ']]
            
        eq_num = 1
        for n_id in sorted(self.nodes.keys()):
            for dof in range(3):
                if self.equation_map[n_id][dof] == 0:  
                    self.equation_map[n_id][dof] = eq_num
                    eq_num += 1
                else:  
                    self.equation_map[n_id][dof] = 0
        self.num_equations = eq_num - 1

    def assemble_matrices(self):
        """Builds local matrices, transforms them, and assembles the Sparse Global Matrix."""
        self.K_global_struct = SparseMatrix(self.num_equations)
        
        for mem_id, mem in self.members.items():
            mat = self.materials[mem['material_id']]
            E, A, I = mat['E'], mat['A'], mat['I']
            
            start_n = self.nodes[mem['start_node']]
            end_n = self.nodes[mem['end_node']]
            
            dx = end_n['X'] - start_n['X']
            dy = end_n['Y'] - start_n['Y']
            
            L = (dx**2 + dy**2)**0.5
            c, s = dx / L, dy / L
            
            # Local Stiffness
            k_loc = Matrix(6, 6)
            k_loc.set_val(0, 0, E * A / L); k_loc.set_val(3, 3, E * A / L)
            k_loc.set_val(0, 3, -E * A / L); k_loc.set_val(3, 0, -E * A / L)
            
            k_loc.set_val(1, 1, 12 * E * I / (L**3)); k_loc.set_val(4, 4, 12 * E * I / (L**3))
            k_loc.set_val(1, 4, -12 * E * I / (L**3)); k_loc.set_val(4, 1, -12 * E * I / (L**3))
            
            k_loc.set_val(1, 2, 6 * E * I / (L**2)); k_loc.set_val(2, 1, 6 * E * I / (L**2))
            k_loc.set_val(1, 5, 6 * E * I / (L**2)); k_loc.set_val(5, 1, 6 * E * I / (L**2))
            k_loc.set_val(4, 2, -6 * E * I / (L**2)); k_loc.set_val(2, 4, -6 * E * I / (L**2))
            k_loc.set_val(4, 5, -6 * E * I / (L**2)); k_loc.set_val(5, 4, -6 * E * I / (L**2))
            
            k_loc.set_val(2, 2, 4 * E * I / L); k_loc.set_val(5, 5, 4 * E * I / L)
            k_loc.set_val(2, 5, 2 * E * I / L); k_loc.set_val(5, 2, 2 * E * I / L)
            
            self.k_local[mem_id] = k_loc
            
            # Rotation Matrix
            R = Matrix(6, 6)
            R.set_val(0, 0, c); R.set_val(1, 1, c); R.set_val(3, 3, c); R.set_val(4, 4, c)
            R.set_val(0, 1, s); R.set_val(3, 4, s)
            R.set_val(1, 0, -s); R.set_val(4, 3, -s)
            R.set_val(2, 2, 1); R.set_val(5, 5, 1)
            
            self.rot_matrices[mem_id] = R
            
            # Global Element Stiffness k_glob = R^T * k_loc * R
            R_T = R.transpose()
            k_glob = R_T.multiply(k_loc).multiply(R)
            self.k_global_element[mem_id] = k_glob
            
            # Assemble into Sparse Global Structure Matrix
            sn_id, en_id = mem['start_node'], mem['end_node']
            G = self.equation_map[sn_id] + self.equation_map[en_id]
            
            for p in range(6):
                for q in range(6):
                    P, Q = G[p], G[q]
                    if P != 0 and Q != 0:  
                        self.K_global_struct.add_val(P - 1, Q - 1, k_glob.get_val(p, q))

    def solve_system(self):
        """Builds load vector and solves for displacements."""
        self.load_vector = [0.0] * self.num_equations
        
        for load in self.loads.values():
            n_id = load['node_id']
            dofs = self.equation_map[n_id] 
            forces = [load['Fx'], load['Fy'], load['Mz']]
            
            for q in range(3):
                Q = dofs[q]
                if Q != 0:
                    self.load_vector[Q - 1] += forces[q]

        self.displacements = self.K_global_struct.solve(self.load_vector)

    def calculate_internal_forces(self):
        """Calculates element end forces in local coordinates."""
        for mem_id, mem in self.members.items():
            sn_id, en_id = mem['start_node'], mem['end_node']
            G = self.equation_map[sn_id] + self.equation_map[en_id]
            
            d_global = [0.0] * 6
            for i in range(6):
                if G[i] != 0:
                    d_global[i] = self.displacements[G[i] - 1] 
                    
            R = self.rot_matrices[mem_id]
            k_loc = self.k_local[mem_id]
            
            d_local = R.multiply(d_global)
            f_local = k_loc.multiply(d_local)
            self.member_forces[mem_id] = f_local