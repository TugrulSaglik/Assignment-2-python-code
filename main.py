"""
Module: main
Purpose: Entry point. Handles file I/O, runs the structural analysis pipeline, 
         and generates standard formatted text reports.
"""

from frame_solver import FrameModel

def parse_input_file(filepath, model):
    """Reads structural data and populates the FrameModel."""
    current_section = None
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'): continue
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].lower()
                continue
            
            parts = [p.strip() for p in line.split('/')]
            if current_section == 'materials':
                model.materials[int(parts[0])] = {'A': float(parts[1]), 'I': float(parts[2]), 'E': float(parts[3])}
            elif current_section == 'nodes':
                model.nodes[int(parts[0])] = {'X': float(parts[1]), 'Y': float(parts[2])}
            elif current_section == 'members':
                model.members[int(parts[0])] = {'start_node': int(parts[1]), 'end_node': int(parts[2]), 'material_id': int(parts[3])}
            elif current_section == 'supports':
                model.supports[int(parts[0])] = {'node_id': int(parts[1]), 'RX': int(parts[2]), 'RY': int(parts[3]), 'RZ': int(parts[4])}
            elif current_section == 'loads':
                model.loads[int(parts[0])] = {'node_id': int(parts[1]), 'Fx': float(parts[2]), 'Fy': float(parts[3]), 'Mz': float(parts[4])}

def write_text_report(model, out_filepath):
    """Writes the comprehensive analysis results to a simple .txt file."""
    with open(out_filepath, 'w') as f:
        f.write("=== STRUCTURAL ANALYSIS REPORT ===\n\n")
        
        f.write("--- EQUATION NUMBERING ---\n")
        for n_id in sorted(model.equation_map.keys()):
            dofs = model.equation_map[n_id]
            f.write(f"Node {n_id:2} | DX: {dofs[0]:2} | DY: {dofs[1]:2} | RZ: {dofs[2]:2}\n")
            
        f.write(f"\nTotal Equations: {model.num_equations}\n\n")
        
        f.write("--- STRUCTURAL DISPLACEMENTS ---\n")
        for i, d in enumerate(model.displacements):
            f.write(f"Eq {i+1:2}: {d:>12.4e}\n")
            
        f.write("\n--- MEMBER END FORCES ---\n")
        for mem_id, forces in model.member_forces.items():
            f.write(f"\nMember {mem_id}:\n")
            f.write(f"  Start -> Axial: {forces[0]:8.2f}, Shear: {forces[1]:8.2f}, Moment: {forces[2]:8.2f}\n")
            f.write(f"  End   -> Axial: {forces[3]:8.2f}, Shear: {forces[4]:8.2f}, Moment: {forces[5]:8.2f}\n")

if __name__ == "__main__":
    model = FrameModel()
    
    # 1. Read Inputs
    parse_input_file('input.txt', model)
    
    # 2. Execute Analysis Pipeline
    model.process_equations()
    model.assemble_matrices()
    model.solve_system()
    model.calculate_internal_forces()
    
    # 3. Generate Output File
    write_text_report(model, 'output_report.txt')
    print("Analysis complete! Check 'output_report.txt' for the results.")