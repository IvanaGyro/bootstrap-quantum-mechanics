# Bootstrap Quantum Mechanics

## Requirements

This project use [Sagemath](https://www.sagemath.org/). Please following the [official docuement to install Sagemath](https://doc.sagemath.org/html/en/installation/index.html)

The code is only tested with Python 3.10. It does not garantee to work with other versions.

### For Windows

If you run Sagemath in Windows Subsystem for Linux (WSL), you may need X11 and VcXsrv to show the output images. Please refer to [this guide](https://www.guide2wsl.com/x11/) to set up the environment.

## Execution

1. Modifiy the maximum depth in `main.sage`:

```python
if __name__ == '__main__':
    find_possible_energy((2, 80), (-0.7, 0.5)) # Maximum depth is 80 by default.
```

If you don't want to wait too long, 20 is a good choice.

```python
if __name__ == '__main__':
    find_possible_energy((2, 20), (-0.7, 0.5))
```

2. Run `sage main.sage`. The program will automatically create `dumps/` and `images/` folders under the root folder of this project.
