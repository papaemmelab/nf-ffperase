import os
import subprocess
import time
from openai import OpenAI

from rich.console import Console
from rich.live import Live
from rich.panel import Panel

console = Console()

dockerfile_path = "Dockerfile.test"

def build_image():
    try:
        # Stream output from docker build command
        proc = subprocess.Popen(
            ["docker", "build", "."],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        output_lines = []
        # Print output line by line in real time
        for line in proc.stdout:
            output_lines.append(line)
            console.print(line, style="cyan")
        proc.wait()
        combined_output = "".join(output_lines)
        return proc.returncode, combined_output
    except Exception as e:
        return 1, str(e)

def ask_llm_to_fix(error_log, dockerfile_content):
    api_key = os.getenv("OPENAI_API_KEY")
    client = OpenAI(api_key=api_key)
    prompt = f"""
        You're a DevOps assistant. The user tried to build a Dockerfile and got this error:
        \n\n
        {error_log}
        \n\n
        Here is the Dockerfile:\n\n{dockerfile_content}\n\nSuggest an improved version of the Dockerfile to fix the error. 
        Please only respond with a valid Dockerfile.
        I need to have installed matplotlib, pandas and pycirclize with any version of python >3.8
    """
    response = client.chat.completions.create(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}]
    )
    return response.choices[0].message.content

# Main loop with visualization of current Dockerfile and output log
for iteration in range(5):
    with open(dockerfile_path, "r") as f:
        dockerfile = f.read()

    # Use Rich Live to update a panel with the Dockerfile content
    with Live(Panel(dockerfile, title="Current Dockerfile", border_style="green"), refresh_per_second=2) as live:
        code, output = build_image()
        live.update(Panel(dockerfile, title="Current Dockerfile", border_style="green"))
        console.print(f"[bold]Iteration {iteration+1}: Building image...[/bold]")
        if code == 0:
            console.print("[bold green]✅ Build successful![/bold green]")
            break
        else:
            console.print("[bold red]❌ Build failed. Asking LLM for help...[/bold red]")
            console.print(Panel(output, title="Error Log", border_style="red"))
            suggestion = ask_llm_to_fix(output, dockerfile)
            console.print(Panel(suggestion, title="LLM Suggestion", border_style="blue"))
            console.print("🔁 Updating Dockerfile...\n")
            with open(dockerfile_path, "w") as f:
                f.write(suggestion)
            # Give a moment to inspect before re-building
            time.sleep(2)