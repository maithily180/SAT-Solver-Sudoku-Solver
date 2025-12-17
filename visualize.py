import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

# Read CSV files with proper dtypes
try:
    bt_data = pd.read_csv('results/backtracking_results.csv')
    sat_data = pd.read_csv('results/sat_results.csv')
    
    # Convert numeric columns from strings to numbers
    bt_data['Time_s'] = pd.to_numeric(bt_data['Time_s'], errors='coerce')
    bt_data['Decisions'] = pd.to_numeric(bt_data['Decisions'], errors='coerce')
    bt_data['Backtracks'] = pd.to_numeric(bt_data['Backtracks'], errors='coerce')
    
    # Handle TIMEOUT in SAT data
    sat_data['Time_s'] = sat_data['Time_s'].replace('TIMEOUT', '30.0')
    sat_data['Time_s'] = pd.to_numeric(sat_data['Time_s'], errors='coerce')
    sat_data['Decisions'] = pd.to_numeric(sat_data['Decisions'], errors='coerce')
    sat_data['UnitProps'] = pd.to_numeric(sat_data['UnitProps'], errors='coerce')
    sat_data['Backtracks'] = pd.to_numeric(sat_data['Backtracks'], errors='coerce')
    sat_data['Backjumps'] = pd.to_numeric(sat_data['Backjumps'], errors='coerce')
    
    print("✓ Data loaded successfully!")
    print(f"  Backtracking: {len(bt_data)} records")
    print(f"  SAT: {len(sat_data)} records")
except FileNotFoundError as e:
    print(f"Error: Could not find CSV files. Make sure you've run the solver first!")
    exit(1)

# Create output directory
import os
os.makedirs('analysis', exist_ok=True)

# ============================================================================
# GRAPH 1: Backtracking Performance
# ============================================================================
def plot_backtracking_performance():
    print("\nGenerating Graph 1: Backtracking Performance...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    bt_grouped = bt_data.groupby('Configuration').agg({
        'Time_s': 'mean',
        'Decisions': 'mean',
        'Backtracks': 'mean'
    }).reset_index().sort_values('Time_s')
    
    colors = ['#2ecc71' if 'MRV' in c else '#e74c3c' for c in bt_grouped['Configuration']]
    
    ax1.barh(bt_grouped['Configuration'], bt_grouped['Time_s'], color=colors, alpha=0.8)
    ax1.set_xlabel('Time (seconds)', fontweight='bold')
    ax1.set_title('Backtracking: Time Comparison', fontweight='bold', fontsize=14)
    ax1.grid(axis='x', alpha=0.3)
    
    for i, (idx, row) in enumerate(bt_grouped.iterrows()):
        ax1.text(row['Time_s'], i, f" {row['Time_s']:.4f}s", va='center', fontsize=9)
    
    ax2.barh(bt_grouped['Configuration'], bt_grouped['Decisions'], color=colors, alpha=0.8)
    ax2.set_xlabel('Number of Decisions', fontweight='bold')
    ax2.set_title('Backtracking: Search Effort', fontweight='bold', fontsize=14)
    ax2.grid(axis='x', alpha=0.3)
    
    for i, (idx, row) in enumerate(bt_grouped.iterrows()):
        ax2.text(row['Decisions'], i, f" {int(row['Decisions'])}", va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('analysis/backtracking_performance.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved: analysis/backtracking_performance.png")
    plt.close()

# ============================================================================
# GRAPH 2: SAT Performance
# ============================================================================
def plot_sat_performance():
    print("\nGenerating Graph 2: SAT Solver Performance...")
    
    sat_success = sat_data[sat_data['Timeout'] == 'NO'].copy()
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    sat_grouped = sat_success.groupby('Configuration')['Time_s'].mean().sort_values()
    colors = ['#27ae60' if i < len(sat_grouped)//2 else '#f39c12' for i in range(len(sat_grouped))]
    
    ax1.barh(sat_grouped.index, sat_grouped.values, color=colors, alpha=0.8)
    ax1.set_xlabel('Time (seconds)', fontweight='bold')
    ax1.set_title('SAT: Time Comparison (Successful)', fontweight='bold', fontsize=13)
    ax1.grid(axis='x', alpha=0.3)
    
    for i, (config, time) in enumerate(sat_grouped.items()):
        ax1.text(time, i, f" {time:.4f}s", va='center', fontsize=9)
    
    decisions = sat_success.groupby('Configuration')['Decisions'].mean().sort_values()
    ax2.barh(decisions.index, decisions.values, color=colors, alpha=0.8)
    ax2.set_xlabel('Number of Decisions', fontweight='bold')
    ax2.set_title('SAT: Search Decisions', fontweight='bold', fontsize=13)
    ax2.grid(axis='x', alpha=0.3)
    
    for i, (config, dec) in enumerate(decisions.items()):
        ax2.text(dec, i, f" {int(dec)}", va='center', fontsize=9)
    
    unit_props = sat_success.groupby('Configuration')['UnitProps'].mean().sort_values()
    ax3.barh(unit_props.index, unit_props.values, color='#3498db', alpha=0.8)
    ax3.set_xlabel('Unit Propagations', fontweight='bold')
    ax3.set_title('SAT: Unit Propagation Count', fontweight='bold', fontsize=13)
    ax3.grid(axis='x', alpha=0.3)
    
    for i, (config, props) in enumerate(unit_props.items()):
        ax3.text(props, i, f" {int(props)}", va='center', fontsize=9)
    
    success_count = len(sat_success['Configuration'].unique())
    timeout_count = len(sat_data[sat_data['Timeout'] == 'YES']['Configuration'].unique())
    
    ax4.bar(['Successful', 'Timeout'], [success_count, timeout_count], 
            color=['#27ae60', '#e74c3c'], alpha=0.8, width=0.5)
    ax4.set_ylabel('Number of Configurations', fontweight='bold')
    ax4.set_title('SAT: Success Rate', fontweight='bold', fontsize=13)
    ax4.set_ylim(0, max(success_count, timeout_count) + 1)
    
    for i, (label, count) in enumerate([('Successful', success_count), ('Timeout', timeout_count)]):
        ax4.text(i, count, str(count), ha='center', va='bottom', fontweight='bold', fontsize=14)
    
    plt.tight_layout()
    plt.savefig('analysis/sat_performance.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved: analysis/sat_performance.png")
    plt.close()

# ============================================================================
# GRAPH 3: Backjumping Analysis
# ============================================================================
def plot_backjumping_analysis():
    print("\nGenerating Graph 3: Backjumping Analysis...")
    
    sat_success = sat_data[sat_data['Timeout'] == 'NO'].copy()
    sat_success = sat_success[sat_success['Backtracks'] > 0]
    
    sat_success['BackjumpRate'] = (sat_success['Backjumps'] / sat_success['Backtracks'] * 100)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    grouped = sat_success.groupby('Configuration').agg({
        'Backtracks': 'mean',
        'Backjumps': 'mean'
    }).reset_index()
    
    x = np.arange(len(grouped))
    width = 0.35
    
    ax1.bar(x - width/2, grouped['Backtracks'], width, label='Total Backtracks', 
            color='#e74c3c', alpha=0.8)
    ax1.bar(x + width/2, grouped['Backjumps'], width, label='Backjumps', 
            color='#27ae60', alpha=0.8)
    
    ax1.set_xlabel('Configuration', fontweight='bold')
    ax1.set_ylabel('Count', fontweight='bold')
    ax1.set_title('Backjumping: Chronological vs Non-Chronological', fontweight='bold', fontsize=13)
    ax1.set_xticks(x)
    ax1.set_xticklabels(grouped['Configuration'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    jump_rate = sat_success.groupby('Configuration')['BackjumpRate'].mean().sort_values()
    
    colors = ['#27ae60' if rate > 50 else '#f39c12' for rate in jump_rate.values]
    ax2.barh(jump_rate.index, jump_rate.values, color=colors, alpha=0.8)
    ax2.set_xlabel('Backjump Rate (%)', fontweight='bold')
    ax2.set_title('% of Non-Chronological Backtracks', fontweight='bold', fontsize=13)
    ax2.grid(axis='x', alpha=0.3)
    ax2.axvline(50, color='red', linestyle='--', alpha=0.5, label='50% threshold')
    ax2.legend()
    
    for i, (config, rate) in enumerate(jump_rate.items()):
        ax2.text(rate, i, f" {rate:.1f}%", va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('analysis/backjumping_analysis.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved: analysis/backjumping_analysis.png")
    plt.close()

# ============================================================================
# GRAPH 4: Unit Propagation Impact
# ============================================================================
def plot_unit_prop_impact():
    print("\nGenerating Graph 4: Unit Propagation Impact...")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    sat_with_up = sat_data[sat_data['Configuration'].str.contains('Unit Prop|All')]
    with_up_time = sat_with_up[sat_with_up['Timeout'] == 'NO']['Time_s'].mean()
    without_up_time = 30.0
    
    categories = ['With Unit\nPropagation', 'Without Unit\nPropagation']
    times = [with_up_time, without_up_time]
    colors = ['#27ae60', '#e74c3c']
    
    bars = ax.bar(categories, times, color=colors, alpha=0.8, width=0.5)
    ax.set_ylabel('Average Time (seconds)', fontweight='bold')
    ax.set_title('Critical Impact of Unit Propagation', fontweight='bold', fontsize=14)
    ax.set_ylim(0, 35)
    
    ax.text(0, with_up_time, f'{with_up_time:.4f}s\n✓ Success', 
            ha='center', va='bottom', fontweight='bold', fontsize=12)
    ax.text(1, without_up_time, f'>30s\n✗ Timeout', 
            ha='center', va='bottom', fontweight='bold', fontsize=12, color='red')
    
    speedup = without_up_time / with_up_time
    ax.annotate(f'{speedup:.0f}x\nSpeedup!', 
                xy=(0.5, (with_up_time + without_up_time)/2),
                xytext=(0.5, 25),
                ha='center',
                fontsize=16,
                fontweight='bold',
                color='#2c3e50',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                arrowprops=dict(arrowstyle='->', lw=2, color='#2c3e50'))
    
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('analysis/unit_propagation_impact.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved: analysis/unit_propagation_impact.png")
    plt.close()

# ============================================================================
# GRAPH 5: Solver Comparison
# ============================================================================
def plot_solver_comparison():
    print("\nGenerating Graph 5: Solver Comparison...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    bt_best = bt_data.loc[bt_data['Time_s'].idxmin()]
    sat_best = sat_data[sat_data['Timeout'] == 'NO'].loc[
        sat_data[sat_data['Timeout'] == 'NO']['Time_s'].idxmin()
    ]
    
    solvers = ['Backtracking\n(Best)', 'SAT\n(Best)']
    times = [bt_best['Time_s'], sat_best['Time_s']]
    colors = ['#3498db', '#9b59b6']
    
    bars = ax1.bar(solvers, times, color=colors, alpha=0.8, width=0.5)
    ax1.set_ylabel('Time (seconds)', fontweight='bold')
    ax1.set_title('Best Configuration: Time Comparison', fontweight='bold', fontsize=13)
    
    for i, (solver, time) in enumerate(zip(solvers, times)):
        ax1.text(i, time, f'{time:.4f}s', ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    ax1.grid(axis='y', alpha=0.3)
    
    decisions = [bt_best['Decisions'], sat_best['Decisions']]
    
    bars = ax2.bar(solvers, decisions, color=colors, alpha=0.8, width=0.5)
    ax2.set_ylabel('Number of Decisions', fontweight='bold')
    ax2.set_title('Best Configuration: Search Effort', fontweight='bold', fontsize=13)
    
    for i, (solver, dec) in enumerate(zip(solvers, decisions)):
        ax2.text(i, dec, f'{int(dec)}', ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('analysis/solver_comparison.png', dpi=300, bbox_inches='tight')
    print("  ✓ Saved: analysis/solver_comparison.png")
    plt.close()

# ============================================================================
# GENERATE ALL GRAPHS
# ============================================================================
print("=" * 70)
print("GENERATING VISUALIZATION AND ANALYSIS")
print("=" * 70)

plot_backtracking_performance()
plot_sat_performance()
plot_backjumping_analysis()
plot_unit_prop_impact()
plot_solver_comparison()

print("\n" + "=" * 70)
print("✓ ALL GRAPHS GENERATED SUCCESSFULLY!")
print("=" * 70)
print(f"\nGraphs saved in: analysis/")
print("Files created:")
print("  1. backtracking_performance.png")
print("  2. sat_performance.png")
print("  3. backjumping_analysis.png")
print("  4. unit_propagation_impact.png")
print("  5. solver_comparison.png")
